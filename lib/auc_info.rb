require 'models'
require 'sequence_dataset'
require 'median'

class AUCInfo
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

  def median_auc(model_name)
    median(aucs_for_model(model_name))
  end

  def median_aucs_by_model
    model_names.map{|model_name|
      [model_name, median_auc(model_name)]
    }.to_h
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
