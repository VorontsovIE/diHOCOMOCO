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

def aucs_for_hocomoco(hocomoco_models, auc_infos_for_uniprot, auc_infos_for_uniprot_not_filtered, quality)
  hocomoco_models.select{|model|
    Models.hocomoco_qualities[model.model_name] == quality
  }.map{|model|
    uniprot = model.uniprot
    good_datasets = auc_infos_for_uniprot[uniprot].datasets
    good_datasets.map{|dataset| # take only good datasets
      {
        dataset: dataset,
        model: model,
        auc: auc_infos_for_uniprot_not_filtered[uniprot].auc_by_model_and_dataset[model][dataset],
        dataset_quality: auc_infos_for_uniprot[uniprot].dataset_qualities[dataset],
      }
    }
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
  auc_infos_for_uniprot_not_filtered = AUCs.load_auc_infos_for_uniprot
  auc_infos_for_uniprot = AUCs.load_auc_infos_for_uniprot(min_weight_for_dataset: 0.65, min_auc_for_model: 0.65)
  # auc_infos_for_uniprot = auc_by_model_and_control('occurences/auc/')

  hocomoco_models = Models.all_models.select{|model|
    model.collection_short_name == 'HL'
  }.select{|model|
    model.species == 'HUMAN'
  }.select{|model|
    auc_infos_for_uniprot[model.uniprot]
  }.reject{|model|
    auc_infos_for_uniprot[model.uniprot].empty? # There are TF datasets but model not necessarily weighted
  }

  $stderr.puts "Total #{hocomoco_models.size} hocomoco models with chipseqs"

  thresholds = (0.0..1.0).step(0.04)

  result_total = []
  ['A', 'B', 'C', 'D'].each do |quality|
    auc_infos = aucs_for_hocomoco(hocomoco_models, auc_infos_for_uniprot, auc_infos_for_uniprot_not_filtered, quality)
    max_vals = auc_infos.map{|model_auc_infos|
      model_auc_infos.map{|h| h[:auc] }.max
    }
    min_vals = auc_infos.map{|model_auc_infos|
      model_auc_infos.map{|h| h[:auc] }.min
    }
    median_vals = auc_infos.map{|model_auc_infos|
      median(model_auc_infos.map{|h| h[:auc] })
    }
    mean_vals = auc_infos.map{|model_auc_infos|
      weighted_sum = model_auc_infos.map{|h| h[:auc] * h[:dataset_quality] }.inject(0.0, &:+)
      total_weight_sum = model_auc_infos.map{|h| h[:dataset_quality] }.inject(0.0, &:+)
      weighted_sum / total_weight_sum
    }

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

desc 'Output information about AUCs'
task :hocomoco_model_AUCs do
  auc_infos_for_uniprot_not_filtered = AUCs.load_auc_infos_for_uniprot
  auc_infos_for_uniprot = AUCs.load_auc_infos_for_uniprot(min_weight_for_dataset: 0.65, min_auc_for_model: 0.65)

  hocomoco_models = Models.all_models.select{|model|
    model.collection_short_name == 'HL'
  }.select{|model|
    model.species == 'HUMAN'
  }.select{|model|
    auc_infos_for_uniprot[model.uniprot]
  }.reject{|model|
    auc_infos_for_uniprot[model.uniprot].empty? # There are TF datasets but model not necessarily weighted
  }.sort

  $stderr.puts "Total #{hocomoco_models.size} hocomoco models with chipseqs"

  puts ['Model', 'Quality', 'Max AUC', 'Min AUC', 'Weighted AUC', 'Weighted AUC (original)'].join("\t")
  hocomoco_models.each do |model|
    uniprot = model.uniprot
    good_datasets = auc_infos_for_uniprot[uniprot].datasets
    auc_infos = good_datasets.map{|dataset|
      {
        dataset: dataset,
        model: model,
        auc: auc_infos_for_uniprot_not_filtered[uniprot].auc_by_model_and_dataset[model][dataset],
        dataset_quality: auc_infos_for_uniprot[uniprot].dataset_qualities[dataset],
      }
    }

    weighted_sum = auc_infos.map{|h| h[:auc] * h[:dataset_quality] }.inject(0.0, &:+)
    total_weight_sum = auc_infos.map{|h| h[:dataset_quality] }.inject(0.0, &:+)
    mean_auc = weighted_sum / total_weight_sum
    max_auc = auc_infos.map{|h| h[:auc] }.max
    min_auc = auc_infos.map{|h| h[:auc] }.min
    auc = auc_infos_for_uniprot[uniprot].weighted_model_aucs[model]

    puts [model.model_name, Models.hocomoco_qualities[model.model_name], max_auc, min_auc, mean_auc, auc].join("\t")
  end
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
