require 'models'
require 'best_models'

task :dataset_qualities do
  auc_infos_for_uniprot = AUCs.load_auc_infos_for_uniprot(min_weight_for_dataset: 0.65,
                                                          min_auc_for_model: 0.65)

  # p SequenceDataset.each_dataset.map(&:name) - auc_infos_for_uniprot.values.flat_map(&:all_datasets)
  # exit
  headers = ['Dataset', 'UniprotID', 'Dataset weight', 'Status']
  puts headers.join("\t")

  SequenceDataset.each_dataset{|dataset|
    uniprot = dataset.uniprot
    dataset_name = dataset.name
    auc_infos = auc_infos_for_uniprot[uniprot]

    if auc_infos.models.empty?
      infos = [dataset_name, uniprot, nil, auc_infos.all_models.empty? ? 'No models passed curration' : 'No models passed wAUC threshold']
    else
      good_dataset = auc_infos.datasets.include?(dataset_name)
      infos = [dataset_name, uniprot, auc_infos.dataset_quality(dataset_name), good_dataset ? 'OK' : 'rejected']
    end

    puts infos.join("\t")
  }
end
