require 'dataset_qualities'
require 'models'
require 'sequence_dataset'
require 'auc_info'
require 'auc_infos_filtering'

desc 'Filter too bad datasets and models'
task :filter_motifs_and_datasets do
  dataset_quality_by_name = DatasetQuality.each_in_xlsx('check_result.xlsx').map{|infos|
    [infos.dataset_name, infos]
  }.to_h

  filtering = AUCInfosFiltering.new(AUCInfo.load_all_infos)
  filtering.remove_bad_datasets_and_models!(dataset_quality_by_name, 0.6)

  collection_perfomances = filtering.model_perfomances_collections_grouped

  collections = filtering.model_names.map{|model_name| Model.get_collection_short_name(model_name) }.uniq.sort

  File.open('collection_perfomances.tsv', 'w') do |fw|

    # Median AUC of best mononucleotide model, Best mononucleotide collection, Best mononucleotide model
    # Median AUC of best dinucleotide model, Best dinucleotide collection, Best dinucleotide model
    headers = ['Uniprot', 'Species', 'Median AUC of best model', 'Best collection', 'Best model', *collections]
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
      best_model, best_auc = best_model_with_auc_by_collection[best_collection]
      fw.puts [uniprot, uniprot[/_(?<species>HUMAN|MOUSE)$/,:species], best_auc, best_collection, best_model, *aucs].join("\t")
    }
  end
end
