require 'models'
require 'sequence_dataset'
require 'median'

class AUCInfo
  attr_reader :auc_by_dataset_and_model, :auc_by_model_and_dataset
  def initialize
    @auc_by_model_and_dataset = Hash.new{|h,k| h[k] = {} }
    @auc_by_dataset_and_model = Hash.new{|h,k| h[k] = {} }
  end

  def dataset_names; @auc_by_dataset_and_model.keys; end
  def model_names; @auc_by_model_and_dataset.keys; end

  def dataset_names_for_model(model_name); @auc_by_model_and_dataset[model_name].keys; end
  def model_names_for_dataset(dataset_name); @auc_by_dataset_and_model[dataset_name].keys; end

  def aucs_for_dataset(dataset_name); @auc_by_dataset_and_model[dataset_name].values; end
  def aucs_for_model(model_name); @auc_by_model_and_dataset[model_name].values; end

  def auc(model_name, dataset_name)
    @auc_by_model_and_dataset[model_name][dataset_name]
  end

  # median AUC over datasets
  def median_auc(model_name)
    median(aucs_for_model(model_name))
  end

  # max AUC over datasets
  def max_auc(model_name)
    aucs_for_model(model_name).max
  end

  # {model => auc}
  def median_aucs_by_model
    model_names.map{|model_name|
      [model_name, median_auc(model_name)]
    }.to_h
  end

  # {model => auc}
  def max_aucs_by_model
    model_names.map{|model_name|
      [model_name, max_auc(model_name)]
    }.to_h
  end

  # {uniprot => {model => auc} }
  def model_perfomances_by_uniprot(aucs_by_model)
    aucs_by_model.group_by{|model_name, auc|
      Model.get_uniprot(model_name)
    }.map{|uniprot, model_auc_pairs|
      [uniprot, model_auc_pairs.to_h]
    }.to_h
  end

  def max_model_perfomances_by_uniprot
    model_perfomances_by_uniprot(max_aucs_by_model)
  end

  def median_model_perfomances_by_uniprot
    model_perfomances_by_uniprot(median_aucs_by_model)
  end

  # {uniprot => {collection => [model, auc] } }
  def model_perfomances_collections_grouped(aucs_by_model)
    model_perfomances_by_uniprot(aucs_by_model).map{|uniprot, auc_by_model|
      best_model_with_auc_by_collection = auc_by_model.group_by{|model_name, auc|
        Model.get_collection_short_name(model_name)
      }.map{|short_collection_id, model_aucs|
        [short_collection_id, model_aucs.max_by{|model_name, auc| auc }]
      }.to_h
      [uniprot, best_model_with_auc_by_collection]
    }.to_h
  end

  def max_model_perfomances_collections_grouped
    model_perfomances_collections_grouped(max_aucs_by_model)
  end

  def median_model_perfomances_collections_grouped
    model_perfomances_collections_grouped(median_aucs_by_model)
  end

  # Dataset removal operations
  def remove_dataset!(dataset_name)
    @auc_by_dataset_and_model.delete(dataset_name)
    @auc_by_model_and_dataset.each_value{|auc_by_dataset|
      auc_by_dataset.delete(dataset_name)
    }
  end

  def remove_datasets!(dataset_names_to_remove)
    dataset_names_to_remove.each{|dataset_name| remove_dataset!(dataset_name) }
  end

  def remove_datasets_wo_models!
    @auc_by_dataset_and_model.delete_if{|dataset_name, auc_by_model|
      auc_by_model.empty?
    }
  end

  # Model removal operations
  def remove_model!(model_name)
    @auc_by_model_and_dataset.delete(model_name)
    @auc_by_dataset_and_model.each_value{|auc_by_model|
      auc_by_model.delete(model_name)
    }
  end

  def remove_models!(models_to_remove)
    models_to_remove.each{|model_name| remove_model!(model_name) }
  end

  def remove_models_wo_datasets!
    @auc_by_model_and_dataset.delete_if{|model_name, auc_by_dataset|
      auc_by_dataset.empty?
    }
  end

  # Information addition
  def add_data_from_file!(auc_filename, model_fullname)
    File.readlines(auc_filename).each{|line|
      dataset_name, auc = line.chomp.split("\t")
      auc = auc.to_f
      @auc_by_model_and_dataset[model_fullname][dataset_name] = auc
      @auc_by_dataset_and_model[dataset_name][model_fullname] = auc
    }
  end

  def self.load_all_infos
    result = AUCInfo.new
    SequenceDataset.each_uniprot do |uniprot|
      Models.all_models_by_uniprot(uniprot).each{|model|
        model_auc_fn = File.join('occurences/auc/', uniprot, "#{model.full_name}.txt")
        result.add_data_from_file!(model_auc_fn, model.full_name)
      }
    end
    result
  end
end
