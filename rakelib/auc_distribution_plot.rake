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
  auc_infos_for_uniprot = AUCs.load_auc_infos_for_uniprot(min_weight_for_dataset: 0.65, min_auc_for_model: 0.65)

  hocomoco_models = Models.all_models.select{|model|
    model.collection_short_name == 'HL'
  }.select{|model|
    model.species == 'HUMAN'
  }.sort

  $stderr.puts "Total #{hocomoco_models.size} hocomoco models with chipseqs"

  thresholds = (0.0..1.0).step(0.04)

  result_total = []
  ['A', 'B', 'C', 'D'].each do |quality|
    model_infos = hocomoco_models.select{|model|
      Models.hocomoco_qualities[model.model_name] == quality
    }.map{|model|
      auc_infos = auc_infos_for_uniprot[model.uniprot]
      aucs = auc_infos.aucs_for_model(model).values
      {
        max_auc: aucs.max,
        min_auc: aucs.min,
        median_auc: aucs.empty? ? nil : median(aucs),
        mean_auc: auc_infos.weighted_auc(model),
      }
    }
    max_vals = model_infos.map{|infos| infos[:max_auc] }
    min_vals =  model_infos.map{|infos| infos[:min_auc] }
    median_vals =  model_infos.map{|infos| infos[:median_auc] }
    mean_vals =  model_infos.map{|infos| infos[:mean_auc] }

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
  auc_infos_for_uniprot = AUCs.load_auc_infos_for_uniprot(min_weight_for_dataset: 0.65, min_auc_for_model: 0.65)

  hocomoco_models = Models.all_models.select{|model|
    model.collection_short_name == 'HL'
  }.select{|model|
    model.species == 'HUMAN'
  }.select{|model|
    auc_infos_for_uniprot[model.uniprot].has_validation?
  }.sort

  $stderr.puts "Total #{hocomoco_models.size} hocomoco models with chipseqs"

  puts ['Model', 'Quality', 'Max AUC', 'Min AUC', 'Weighted AUC', 'Weighted AUC (original)'].join("\t")
  hocomoco_models.each do |model|
    uniprot = model.uniprot
    auc_infos = auc_infos_for_uniprot[uniprot]
    good_datasets = auc_infos.datasets
    aucs = good_datasets.map{|dataset|
      auc_infos_for_uniprot[uniprot].auc(model, dataset)
    }

    max_auc = aucs.max
    min_auc = aucs.min
    auc = auc_infos_for_uniprot[uniprot].weighted_auc(model)

    puts [model.model_name, Models.hocomoco_qualities[model.model_name], max_auc, min_auc, auc, (auc_infos.models.include?(model) ? auc : nil)].join("\t")
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
