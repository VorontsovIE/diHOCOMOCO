require 'dataset_qualities'
require 'models'
require 'sequence_dataset'
require 'auc_info'
require 'auc_infos_filtering'
require 'quality_assessor'

desc 'Filter too bad datasets and models'
task :filter_motifs_and_datasets do
  dataset_quality_by_name = DatasetQuality.each_in_xlsx('check_result.xlsx').map{|infos|
    [infos.dataset_name, infos]
  }.to_h

  filtering = AUCInfosFiltering.new(AUCInfo.load_all_infos)
  filtering.remove_bad_datasets_and_models!(dataset_quality_by_name, 0.6)

  quality_assessor = QualityAssessor.new(filtering)

  collection_perfomances = filtering.model_perfomances_collections_grouped

  collections = filtering.model_names.map{|model_name| Model.get_collection_short_name(model_name) }.uniq.sort
  mono_collections = ['HL', 'HO', 'SR', 'JA', 'SE', 'SMF', 'SMS', 'CM']
  di_collections = ['SDF', 'SDS', 'CD']
  hocomoco_qualities = File.readlines('hocomoco_qualities.tsv').map{|line| line.chomp.split("\t") }.to_h

  File.open('collection_perfomances.tsv', 'w') do |fw|

    headers = ['Uniprot', 'Species',
#              'num_datasets', 'num_survived_models',
              'Median AUC of best model', 'Best collection', 'Best model', 'Quality',
              'Mono or dinucleotide win',
              'Median AUC of best mononucleotide model', 'Best mononucleotide collection', 'Best mononucleotide model', 'Quality of mononucleotide model',
              'Median AUC of best dinucleotide model', 'Best dinucleotide collection', 'Best dinucleotide model', 'Quality of dinucleotide model',
              *collections]
    fw.puts headers.join("\t")
    collection_perfomances.each{|uniprot, best_model_with_auc_by_collection|
      aucs = collections.map{|collection|
        model_name, auc = best_model_with_auc_by_collection[collection]
        auc
      }
      best_collection = collections.max_by{|collection|
        model_name, auc = best_model_with_auc_by_collection[collection]
        auc || 0
      }
      mono_best_collection = mono_collections.max_by{|collection|
        model_name, auc = best_model_with_auc_by_collection[collection]
        auc || 0
      }
      di_best_collection = di_collections.max_by{|collection|
        model_name, auc = best_model_with_auc_by_collection[collection]
        auc || 0
      }
      best_model, best_auc = best_model_with_auc_by_collection[best_collection]
      mono_best_model, mono_best_auc = best_model_with_auc_by_collection[mono_best_collection]
      di_best_model, di_best_auc = best_model_with_auc_by_collection[di_best_collection]

      # ## num_datasets = quality_assessor.num_datasets(best_model)
      # ## num_survived_models = quality_assessor.num_survived_models(best_model)
      # num_datasets_pass_highquality_auc = num_datasets_passing_auc(model_name, 0.9)
      # num_datasets_pass_optimal_auc = num_datasets_passing_auc(model_name, 0.7)
      # num_datasets_pass_minimal_auc = num_datasets_passing_auc(model_name, 0.6)
      quality = quality_assessor.calculate_quality(best_model)

      mono_best_quality = mono_best_model ? quality_assessor.calculate_quality(mono_best_model) : nil
      di_best_quality = di_best_model ? quality_assessor.calculate_quality(di_best_model) : nil


      # add info about number of hocomoco models in test
      fw.puts [uniprot, uniprot[/_(?<species>HUMAN|MOUSE)$/,:species],
 #             num_datasets, num_survived_models,
              best_auc, best_collection, best_model, quality,
              mono_collections.include?(best_collection) ? 'Mono' : 'Di',
              mono_best_auc, mono_best_collection, mono_best_model, mono_best_quality,
              di_best_auc, di_best_collection, di_best_model, di_best_quality,
              *aucs].join("\t")
    }
  end
end
