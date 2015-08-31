require 'auc_by_model_and_control'
require 'median'

def cumulative_distribution(thresholds, data)
  thresholds.map{|t|
    data.count{|x| x <= t }
  }
end

def histogram(thresholds, data)
  thresholds.each_cons(2).map{|t1,t2|
    data.count{|x| (t1...t2).include?(x) }
  }
end

def load_hocomoco_qualities(filename)
  File.readlines(filename).drop(1).map{|line|
    line.chomp.split("\t")
  }.map{|row|
    [row[0],row[2]]
  }.to_h
end

def aucs_for_hocomoco(hocomoco_models, aucs_by_model_and_control, hocomoco_qualities, quality)
  hocomoco_models.select{|model|
    hocomoco_qualities[model.model_name] == quality
  }.map{|model|
    aucs_by_model_and_control[model.full_name].values
  }
end

task :auc_distribution_plot do
  auc_by_control_and_model = auc_by_control_and_model('occurences/auc/')
  max_aucs_for_uniprots = SequenceDataset.each_uniprot.map{|uniprot|
    SequenceDataset.each_for_uniprot(uniprot).map{|control|
      auc_by_control_and_model[control.name].map{|model, auc| auc }.max
    }.max
  }

  minimax_aucs_for_uniprots = SequenceDataset.each_uniprot.map{|uniprot|
    SequenceDataset.each_for_uniprot(uniprot).map{|control|
      auc_by_control_and_model[control.name].map{|model, auc| auc }.max
    }.min
  }

  minimin_aucs_for_uniprots = SequenceDataset.each_uniprot.map{|uniprot|
    SequenceDataset.each_for_uniprot(uniprot).map{|control|
      auc_by_control_and_model[control.name].map{|model, auc| auc }.min
    }.min
  }

  maximin_aucs_for_uniprots = SequenceDataset.each_uniprot.map{|uniprot|
    SequenceDataset.each_for_uniprot(uniprot).map{|control|
      auc_by_control_and_model[control.name].map{|model, auc| auc }.min
    }.max
  }


  thresholds = (0.0..1.0).step(0.01)

  max_auc_distribution = cumulative_distribution(thresholds, max_aucs_for_uniprots)
  minimax_auc_distribution = cumulative_distribution(thresholds, minimax_aucs_for_uniprots)
  maximin_auc_distribution = cumulative_distribution(thresholds, maximin_aucs_for_uniprots)
  minimin_auc_distribution = cumulative_distribution(thresholds, minimin_aucs_for_uniprots)

  headers_1 = ['Threshold', '# of factors which have metric (specified below) over datasets and models less or equal to threshold']
  headers_2 = ['Threshold', 'AUC max over models and datasets', 'AUC max over models, min over datasets', 'AUC min over models, max over datasets', 'AUC min over models, min over datasets']
  result = thresholds.zip(max_auc_distribution, minimax_auc_distribution, maximin_auc_distribution, minimin_auc_distribution).map{|*args|
    args.join("\t")
  }.join("\n")
  File.write 'auc_distribution.tsv', headers_1.join("\t") + "\n" + headers_2.join("\t") + "\n" + result
end

task :hocomoco_auc_distribution do
  aucs_by_model_and_control = auc_by_model_and_control('occurences/auc/')

  hocomoco_models = Models.all_models.select{|m|
    m.collection_short_name == 'HL'
  }.select{|m|
    m.species == 'HUMAN'
  }.reject{|m|
    aucs_by_model_and_control[m.full_name].empty?
  }

  thresholds = (0.0..1.0).step(0.04)
  hocomoco_qualities = load_hocomoco_qualities('hocomoco_genes_infos.csv')

  result_total = []
  ['A', 'B', 'C', 'D'].each do |quality|
    max_vals = aucs_for_hocomoco(hocomoco_models, aucs_by_model_and_control, hocomoco_qualities, quality).map(&:max)
    min_vals = aucs_for_hocomoco(hocomoco_models, aucs_by_model_and_control, hocomoco_qualities, quality).map(&:min)
    mean_vals = aucs_for_hocomoco(hocomoco_models, aucs_by_model_and_control, hocomoco_qualities, quality).map{|vals| mean(vals) }
    median_vals = aucs_for_hocomoco(hocomoco_models, aucs_by_model_and_control, hocomoco_qualities, quality).map{|vals| median(vals) }

    result = []
    result << ["Quality: #{quality}", '', '', '', '']
    result << ['Threshold', '# of models which have metric (specified below) less or equal to threshold', '', '', '']
    result << ['', '', '', '', '']
    result << ['Threshold', 'Max AUC', 'Min AUC', 'Median AUC', 'Mean AUC']
    result += thresholds.drop(1).zip( \
                histogram(thresholds, max_vals), \
                histogram(thresholds, min_vals), \
                histogram(thresholds, median_vals), \
                histogram(thresholds, mean_vals) \
              ).to_a
    result_total += result.transpose
    result_total << [''] * result_total.last.size
  end
  File.write 'hocomoco_AUC_distribution.tsv', result_total.transpose.map{|row| row.join("\t") }.join("\n")
end

task :num_removed_datasets_by_threshold do
  auc_by_control_and_model = auc_by_control_and_model('occurences/auc/')
  max_aucs_for_datasets = SequenceDataset.each_dataset.map{|dataset|
    auc_by_control_and_model[dataset.name].values
  }.map(&:max)
  thresholds = (0.0..1.0).step(0.01)
  max_aucs_distribution = cumulative_distribution(thresholds, max_aucs_for_datasets)
  File.write 'max_aucs.tsv', thresholds.zip(max_aucs_distribution).map{|row| row.join("\t") }.join("\n")
end
