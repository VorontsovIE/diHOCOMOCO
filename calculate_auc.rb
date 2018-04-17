####  Input:
####    >controlName:seqLength
####    pvalue	position	orientation
####    >controlName:seqLength
####    pvalue	position	orientation
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

def correct_pvalues(pvalues, median_length:, model_length:)
  pvalues.map{|pvalue|
    1.0 - (1.0 - pvalue) ** (2 * (median_length - model_length + 1))
  }
end

def tpr_by_fpr(roc_curve, fpr_list)
  points = roc_curve.sort_by{|point| point.fpr }
  fpr_list.map{|fpr|
    points.bsearch{|point| fpr <= point.fpr }.tpr
  }
end

model_length = Integer(ARGV[0])

# prints pvalues or ROC curve instead(!) of AUC
print_pvalues = ARGV.delete('--print-pvalues')
print_corrected_pvalues = ARGV.delete('--print-corrected-pvalues')
print_roc = ARGV.delete('--print-roc')
print_logroc = ARGV.delete('--print-logroc')
if ARGV.include?('--print-tpr-at')
  fpr_values = ARGV.delete_at(ARGV.index('--print-tpr-at') + 1).split(',').map{|pval| Float(pval) }
  ARGV.delete('--print-tpr-at')
  print_tpr = true
  raise  if fpr_values.empty?
end

$stdin.each_line.lazy.map(&:chomp).each_slice(2).chunk{|seqName, hitInfo|
  seqName[1..-1].split(':').first # controlName
}.each{|controlName, controlSequencesIter|
  controlSeqLengths = []
  pvalues = []
  controlSequencesIter.each{|seqName, hitInfo|
    controlSeqLengths << seqName.split(':').last.to_i
    pvalue = hitInfo.split("\t").first.to_f
    pvalues << pvalue
  }

  corrected_pvalues = correct_pvalues(pvalues, median_length: median(controlSeqLengths), model_length: model_length)
  roc = roc_curve(corrected_pvalues)
  logroc_xy = roc.map{|point| [point.fpr <= 0 ? 0 : Math.log(point.fpr), point.tpr] }.drop(1)  # we drop point (log 0; 0) and all points which got negative corrected pvalue due to floating point errors

  if print_pvalues
    puts(pvalues)
    next
  end

  if print_corrected_pvalues
    puts(corrected_pvalues)
    next
  end

  if print_roc
    puts(roc.map{|point| "#{point.fpr}\t#{point.tpr}" })
    next
  end

  if print_logroc
    puts(logroc_xy.map{|logfpr, tpr| "#{logfpr}\t#{tpr}" })
    next
  end

  if print_tpr
    puts(tpr_by_fpr(roc, fpr_values).join("\t"))
    next
  end

  roc_auc = calculate_auc(roc)
  logroc_auc = calculate_auc_by_xy(logroc_xy)

  puts [controlName, roc_auc, logroc_auc].join("\t")

#  auprc = calculate_auc_by_xy(precision_recall_curve(corrected_pvalues, 4**model_length))
#  puts [controlName, auprc].join("\t")
}
