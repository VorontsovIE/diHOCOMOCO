####  Input:
####    >controlName:seqLength
####    score	position	orientation
####    >controlName:seqLength
####    score	position	orientation
####  each control has several sequences listed consequently
#### Example:
## >ANDR_HUMAN.PEAKS033090.pics.control:100
## 4.386433799120242       5       -
## >ANDR_HUMAN.PEAKS033090.pics.control:113
## 6.345827621873722       173     +

$:.unshift File.join(__dir__, 'lib')
require 'median'
require 'sequence_dataset'
require 'rake/ext/string'
require 'roc_curve'

def read_bsearch_table(filename)
  File.readlines(filename).map{|l| l.chomp.split("\t").map(&:to_f) }
end

def pvalue_by_score(requested_score, bsearch_table)
  (bsearch_table.bsearch{|score, pvalue| score >= requested_score } || bsearch_table.last).last
end

def correct_pvalues(pvalues, median_length:, model_length:)
  pvalues.map{|pvalue|
    1.0 - (1.0 - pvalue) ** (2 * (median_length - model_length + 1))
  }
end

model_length = Integer(ARGV[0])
thresholds_fn = ARGV[1]

# prints pvalues or ROC curve instead(!) of AUC
print_pvalues = ARGV.delete('--print-pvalues')
print_corrected_pvalues = ARGV.delete('--print-corrected-pvalues')
print_roc = ARGV.delete('--print-roc')
print_logroc = ARGV.delete('--print-logroc')

threshold_pvalue_table = read_bsearch_table(thresholds_fn)

$stdin.each_line.lazy.map(&:chomp).each_slice(2).chunk{|seqName, hitInfo|
  seqName[1..-1].split(':').first # controlName
}.each{|controlName, controlSequencesIter|
  controlSeqLengths = []
  pvalues = []
  controlSequencesIter.each{|seqName, hitInfo|
    controlSeqLengths << seqName.split(':').last.to_i
    score = hitInfo.split("\t").first.to_f
    pvalue = pvalue_by_score(score, threshold_pvalue_table)
    pvalues << pvalue
  }

  corrected_pvalues = correct_pvalues(pvalues, median_length: median(controlSeqLengths), model_length: model_length)
  roc = roc_curve(corrected_pvalues)
  logroc_xy = roc.map{|point| [point.fpr == 0 ? 0 : Math.log(point.fpr), point.tpr] }.drop(1)  # we drop point (log 0; 0)

  puts(pvalues) and next  if print_pvalues
  puts(corrected_pvalues) and  next  if print_corrected_pvalues
  puts(roc.map{|point| "#{point.fpr}\t#{point.tpr}" }) and next  if print_roc
  puts(logroc_xy.map{|logfpr, tpr| "#{logfpr}\t#{tpr}" }) and next  if print_logroc

  roc_auc = calculate_auc(roc)
  logroc_auc = calculate_auc_by_xy(logroc_xy)
  puts [controlName, roc_auc, logroc_auc].join("\t")

#  auprc = calculate_auc_by_xy(precision_recall_curve(corrected_pvalues, 4**model_length))
#  puts [controlName, auprc].join("\t")
}
