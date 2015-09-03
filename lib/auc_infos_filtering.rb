require 'dataset_qualities'
require 'median'
require 'auc_info'
require 'forwardable'

# Class for filtering aggregated information about AUC for different models and datasets and further aggregation of such information
class AUCInfosFiltering
  attr_reader :auc_infos
  extend Forwardable
  def_delegators :@auc_infos, :remove_models!, :remove_datasets!, :remove_models_wo_datasets!, :remove_datasets_wo_models!
  def_delegators :@auc_infos, :auc
  def_delegators :@auc_infos, :model_names, :dataset_names, :aucs_for_dataset, :aucs_for_model
  def_delegators :@auc_infos, :dataset_names_for_model, :model_names_for_dataset
  def_delegators :@auc_infos, :model_perfomances_by_uniprot, :model_perfomances_collections_grouped
  def_delegators :@auc_infos, :max_model_perfomances_by_uniprot, :max_model_perfomances_collections_grouped
  def_delegators :@auc_infos, :median_model_perfomances_by_uniprot, :median_model_perfomances_collections_grouped

  def initialize(auc_infos)
    @auc_infos = auc_infos
    @dataset_quality_by_name = DatasetQuality.each_in_xlsx('check_result.xlsx').map{|infos|
      [infos.dataset_name, infos]
    }.to_h
  end

  def remove_bad_datasets_and_models!(min_auc)
    # bad datasets are those not passing threshold
    bad_datasets = datasets_not_passing_auc_check(min_auc)
    remove_datasets!(bad_datasets)
    remove_models_wo_datasets!

    # # remove models created by bad datasets
    # models_by_bad_datasets = bad_datasets.flat_map{|dataset_name|
    #   @dataset_quality_by_name[dataset_name].model_names
    # }
    # remove_models!(models_by_bad_datasets)
    # remove_datasets_wo_models!
  end

  def datasets_not_passing_auc_check(min_auc)
    dataset_names.select{|dataset_name|
      aucs = aucs_for_dataset(dataset_name)
      num_models_passing_auc = aucs.count{|auc| auc >= min_auc }
      aucs.empty? || num_models_passing_auc < 1 || (num_models_passing_auc <= 1 && aucs.size >= 2)
    }
  end
end
