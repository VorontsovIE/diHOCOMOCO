require 'median'
require 'models'

# AUCs represents AUC values for a single TF (several models, several datasets; AUC for each model over each dataset)
class AUCs
  attr_reader :auc_by_model_and_dataset

  # {model => {dataset => auc}}
  # model should be an instance of class Model (see #best_model_among_collections)
  def initialize(auc_by_model_and_dataset)
    @auc_by_model_and_dataset = auc_by_model_and_dataset.map{|model, dataset_aucs|
      [model, dataset_aucs.sort_by{|dataset, auc| -auc }.to_h]
    }.reject{|model, dataset_aucs|
      dataset_aucs.empty?
    }.to_h
  end

  def auc_by_dataset_and_model
    @auc_by_dataset_and_model ||= begin
      result = {}
      auc_by_model_and_dataset.each{|model, auc_by_dataset|
        auc_by_dataset.each{|dataset, auc|
          result[dataset] ||= {}
          result[dataset][model] = auc
        }
      }
      result
    end
  end

  def to_s
    "<#{auc_by_model_and_dataset}>"
  end

  def ==(other)
    other.is_a?(self.class) && auc_by_model_and_dataset == other.auc_by_model_and_dataset
  end

  def empty?
    models.empty? || datasets.empty?
  end

  def models
    auc_by_model_and_dataset.keys.sort
  end

  def datasets
    @datasets ||= auc_by_model_and_dataset.values.flat_map(&:keys).uniq.sort
  end

  def dataset_qualities
    @dataset_qualities ||= datasets.map{|dataset|
      [dataset, mean(auc_by_dataset_and_model[dataset].values)]
    }.to_h
  end

  def weighted_model_aucs
    @weighted_model_aucs ||= begin
      quality_norm_factor = dataset_qualities.values.inject(0.0, &:+)
      models.map{|model|
        weighted_auc = datasets.map{|dataset|
          auc_by_model_and_dataset[model][dataset] * dataset_qualities[dataset]
        }.inject(0.0, &:+) / quality_norm_factor
        [model, weighted_auc]
      }.sort_by{|model, auc| -auc }.to_h
    end
  end

  def best_model_among_collections(collections, banned: [])
     best_model, best_auc = weighted_model_aucs.select{|model, auc|
      collections.include?(model.collection_short_name)
    }.reject{|model, auc|
      banned.include?(model)
    }.max_by{|model, auc| auc }
    best_model
  end

  def without_bad_datasets(min_auc)
    datasets_to_retain = datasets.select{|dataset|
      dataset_qualities[dataset] >= min_auc
    }
    auc_by_model_and_dataset_retain = auc_by_model_and_dataset.map{|model, auc_by_dataset|
      auc_by_dataset_retain = datasets_to_retain.map{|dataset|
        [dataset, auc_by_dataset[dataset]]
      }.to_h
      [model, auc_by_dataset_retain]
    }.to_h
    self.class.new(auc_by_model_and_dataset_retain)
  end

  def without_bad_models(min_auc)
    auc_by_model_and_dataset_retain = auc_by_model_and_dataset.select{|model, auc_by_dataset|
      weighted_model_aucs[model] >= min_auc
    }.to_h
    self.class.new(auc_by_model_and_dataset_retain)
  end

  # Loads data for a single TF
  def self.from_folder(glob)
    auc_by_model_and_dataset = FileList[glob].map{|fn|
      model = Model.new_by_name(fn.pathmap('%n'))
      auc_by_dataset = File.readlines(fn).map{|line|
        dataset, auc = line.chomp.split("\t")
        [dataset, auc.to_f]
      }.to_h
      [model, auc_by_dataset]
    }.to_h
    self.new(auc_by_model_and_dataset)
  end

  # Load hash {uniprot => AUCs} with iteratively filtered datasets and models
  def self.load_auc_infos_for_uniprot(min_weight_for_dataset:, min_auc_for_model:)
    uniprots = FileList['occurences/auc/*'].pathmap('%n')

    uniprots.map{|uniprot|
      [uniprot, self.from_folder("occurences/auc/#{uniprot}/*.txt")]
    }.map{|uniprot, auc_infos|
      auc_infos_prev = nil
      while auc_infos != auc_infos_prev
        auc_infos_prev = auc_infos
        auc_infos = auc_infos.without_bad_datasets(min_weight_for_dataset).without_bad_models(min_auc_for_model)
      end
      [uniprot, auc_infos]
    }.reject{|uniprot,auc_infos|
      auc_infos.empty?
    }.to_h
  end
end
