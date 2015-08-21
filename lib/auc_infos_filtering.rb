require 'dataset_qualities'
require 'models'
require 'sequence_dataset'
require 'median'
require 'auc_info'
require 'forwardable'

# Class for filtering aggregated information about AUC for different models and datasets and further aggregation of such information
class AUCInfosFiltering
  attr_reader :auc_infos
  extend Forwardable
  def_delegators :@auc_infos, :remove_models!, :remove_datasets!, :remove_models_wo_datasets!, :remove_datasets_wo_models!
  def_delegators :@auc_infos, :model_names, :dataset_names, :aucs_for_dataset

  def initialize(auc_infos)
    @auc_infos = auc_infos
  end

  def remove_bad_datasets_and_models!(dataset_quality_by_name, min_auc)
    # Models were manually currated, some of models are inquality.
    # One should remove models not passed curration step before dataset assessment
    remove_models!(models_failed_curration(dataset_quality_by_name))
    # remove_datasets_wo_models!

    # bad datasets are those not passing threshold
    bad_datasets = datasets_not_passing_auc_check(min_auc)
    remove_datasets!(bad_datasets)
    remove_models_wo_datasets!

    # remove models created by bad datasets
    models_by_bad_datasets = bad_datasets.flat_map{|dataset_name|
      dataset_quality_by_name[dataset_name].model_names
    }
    remove_models!(models_by_bad_datasets)
    remove_datasets_wo_models!
  end


  def models_failed_curration(dataset_quality_by_name)
    dataset_quality_by_name.each_value.reject(&:pass_quality_control?).flat_map(&:model_names)
  end

  def datasets_not_passing_auc_check(min_auc)
    dataset_names.select{|dataset_name|
      aucs_for_dataset(dataset_name).empty? || aucs_for_dataset(dataset_name).max < min_auc
    }
  end

  # {uniprot => {model => auc} }
  def model_perfomances_by_uniprot
    auc_infos.median_aucs_by_model.group_by{|model_name, auc|
      Model.get_uniprot(model_name)
    }.map{|uniprot, model_auc_pairs|
      [uniprot, model_auc_pairs.to_h]
    }.to_h
  end

  # {uniprot => {collection => [model, auc] } }
  def model_perfomances_collections_grouped
    model_perfomances_by_uniprot.map{|uniprot, auc_by_model|
      best_model_with_auc_by_collection = auc_by_model.group_by{|model_name, auc|
        Model.get_collection_short_name(model_name)
      }.map{|short_collection_id, model_aucs|
        [short_collection_id, model_aucs.max_by{|model_name, auc| auc }]
      }.to_h
      [uniprot, best_model_with_auc_by_collection]
    }.to_h
  end
end