def load_threshold_pvalue_list(fn)
  File.readlines(fn).map{|l| l.chomp.split("\t").map{|x| Float(x) } }
end

def thresholds_by_pvalues_bsearch(threshold_pvalue_list, requested_pvalues)
  requested_pvalues.map{|requested_pvalue|
    ind = threshold_pvalue_list.bsearch_index{|threshold, pvalue| pvalue <= requested_pvalue }
    if threshold_pvalue_list[ind][1] == requested_pvalue
      threshold_pvalue_list[ind][0]
    elsif ind == 0
      threshold_pvalue_list.first[0]
    elsif !ind
      threshold_pvalue_list.last[0]
    else
      (threshold_pvalue_list[ind][0] + threshold_pvalue_list[ind - 1][0]) / 2.0
    end
  }
end

def load_thresholds_by_model(folder, species, arity, requested_pvalues)
  Dir.glob("#{folder}/pwm/*").map{|fn|
    model = File.basename(fn, File.extname(fn))
    threshold_pvalue_list = load_threshold_pvalue_list("#{folder}/thresholds/#{model}.thr")
    threshold_by_pvalue = thresholds_by_pvalues_bsearch(threshold_pvalue_list, requested_pvalues)
    rounded_thresholds = threshold_by_pvalue.map{|pvalue, threshold|
      [pvalue, threshold.round(6)]
    }.to_h
    [model, rounded_thresholds]
  }.to_h
end

def save_standard_thresholds!(filename, infos_for_motifs, requested_pvalues)
  header = ['# P-values', *requested_pvalues]
  matrix = infos_for_motifs.map{|motif_infos|
    [motif_infos[:name], *motif_infos[:standard_thresholds].values_at(*requested_pvalues)]
  }
  File.write(filename, [header, *matrix].map{|row| row.join("\t") }.join("\n"))
end
