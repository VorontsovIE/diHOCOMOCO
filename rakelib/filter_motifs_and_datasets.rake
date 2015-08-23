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
  hocomoco_qualities = File.readlines('hocomoco_genes_infos.csv').map{|line| line.chomp.split("\t") }.map{|row| [row[0],row[2]] }.to_h

  File.open('collection_perfomances.tsv', 'w') do |fw|

    # Median AUC of best mononucleotide model, Best mononucleotide collection, Best mononucleotide model
    # Median AUC of best dinucleotide model, Best dinucleotide collection, Best dinucleotide model
    headers = ['Uniprot', 'Species', 'num_datasets', 'num_survived_models', 'Median AUC of best model', 'Best collection', 'Best model', 'Quality', *collections]
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

      num_datasets = filtering.dataset_names_for_model(best_model).size
      is_hocomoco_model = best_model.match(/~HL~/)
      is_chipseq_model = best_model.match(/~CM~|~CD~/)
      survived_models_by_uniprot = filtering.model_names.group_by{|model_name|
        Model.get_uniprot(model_name)
      }
      num_survived_models = survived_models_by_uniprot[Model.get_uniprot(best_model)].size
      if is_hocomoco_model
        quality = hocomoco_qualities[Model.get_original_model_name(best_model)]
      else
        if num_datasets >= 2
          quality = 'A'
        elsif num_datasets == 1
          if is_chipseq_model
            if num_survived_models > 1
              quality = 'C'
            else
              quality = 'D'
            end
          else
            quality = 'B'
          end
        else
          raise 'Impossible: best model doesn\'t pass filters'
        end
      end

      fw.puts [uniprot, uniprot[/_(?<species>HUMAN|MOUSE)$/,:species], num_datasets, num_survived_models, best_auc, best_collection, best_model, quality, *aucs].join("\t")
    }
  end
end
