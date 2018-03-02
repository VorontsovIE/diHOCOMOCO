require 'median'
require 'models'

def motifs_for_slice(slice_fn)
  banlist = File.readlines('curation/banlist_mono.txt').map(&:strip)
  motif_final = File.basename(slice_fn, '.txt')
  slice_type = motif_final.split('.').last[0]
  semiuniprot = motif_final.split('.').first  # without species part
  all_motifs = File.readlines(slice_fn).map(&:chomp)

  additional_motifs = all_motifs.select{|motif|
    motif.match(/~AD~/)
  }.map{|motif| motif.split('~').last }

  motifs = all_motifs.reject{|motif|
    motif.match(/~AD~/)
  }
  hocomoco_motifs = motifs.select{|motif|
    motif.match(/\.H10MO\./)
  }
  chipseq_motifs = motifs.reject{|motif|
    motif.match(/\.H10MO\./)
  }

  if hocomoco_motifs.empty? && slice_type == 'M'
    hocomoco_motifs = Dir.glob("models/pcm/mono/hocomoco_legacy/#{semiuniprot}_*.H10MO.*.pcm")
                         .map{|fn| File.basename(fn, '.pcm') }
                         .reject{|motif| banlist.include?(motif) }
  end
  { chipseq_motifs: chipseq_motifs, hocomoco_motifs: hocomoco_motifs, additional_motifs: additional_motifs }
end

def motifs_for_di_slice(slice_fn)
  banlist = File.readlines('curation/banlist_di.txt').map(&:strip)
  motif_final = File.basename(slice_fn, '.txt')
  slice_type = motif_final.split('.').last[0]
  semiuniprot = motif_final.split('.').first  # without species part
  all_motifs = File.readlines(slice_fn).map(&:chomp)

  additional_motifs = all_motifs.select{|motif|
    motif.match(/~DIAD~/)
  }.map{|motif| motif.split('~').last }

  motifs = all_motifs.reject{|motif|
    motif.match(/~DIAD~/)
  }
  hocomoco_motifs = motifs.select{|motif|
    motif.match(/\.H10DI\./)
  }
  chipseq_motifs = motifs.reject{|motif|
    motif.match(/\.H10DI\./)
  }.reject{|motif|
    banlist.include?(motif)
  }

  if hocomoco_motifs.empty? && slice_type == 'M'
    hocomoco_motifs = Dir.glob("models/pcm/di/hocomoco_legacy/#{semiuniprot}_*.H10DI.*.dpcm")
                         .map{|fn| File.basename(fn, '.dpcm') }
                         .reject{|motif| banlist.include?(motif) }
  end
  { chipseq_motifs: chipseq_motifs, hocomoco_motifs: hocomoco_motifs, additional_motifs: additional_motifs }
end

# An instance of class AUCs represents AUC values for a single TF
#
# `auc_by_model_and_dataset` is a 2D-hash for AUC for each model-dataset pair
# `models` and `datasets` represent subsets of "good" models and "good" datasets
#   which are treated as reliable enough to be used in mean AUC assessment.
class AUCs
  attr_reader :models, :datasets

  attr_reader :auc_by_model_and_dataset
  protected :auc_by_model_and_dataset

  # {model => {dataset => auc}}
  # model should be an instance of class Model (see #best_model_among_collections)
  #
  # Models and datasets can be specified explicitly to specify good datasets and models
  #   not to use other datasets/models in dataset's weight and model's weighted-AUC assesment
  def initialize(auc_by_model_and_dataset, models: nil, datasets: nil)
    @auc_by_model_and_dataset = auc_by_model_and_dataset

    @models = models || @auc_by_model_and_dataset.keys.sort
    @datasets = datasets || @auc_by_model_and_dataset.values.flat_map(&:keys).uniq.sort
  end

  # All models which were benchmarked, including models which were rejected
  # Use with attention
  def all_models
    @all_models ||= @auc_by_model_and_dataset.keys
  end

  # All models which were benchmarked, including models which were rejected
  # Use with attention
  def all_datasets
    @all_datasets ||= @auc_by_model_and_dataset.values.flat_map(&:keys).uniq
  end

  def to_s
    ["models: #{@models}", "datasets: #{@datasets}", "AUCs: #{@auc_by_model_and_dataset}"].join("\n")
  end
  def inspect; to_s; end

  def ==(other)
    other.is_a?(self.class) && \
    auc_by_model_and_dataset == other.auc_by_model_and_dataset && \
    models == other.models && \
    datasets == other.datasets
  end

  def has_validation?
    !datasets.empty?
  end

  def best_auc(model)
    datasets.map{|ds| auc(model, ds) }.max
  end

  # {dataset => auc}
  def aucs_for_model(model)
    auc_by_model_and_dataset[model].select{|dataset, auc|
      datasets.include?(dataset)
    }
  end

  # {model => auc}
  def aucs_for_dataset(dataset)
    models.map{|model|
      [model, auc_by_model_and_dataset[model][dataset]]
    }.to_h
  end

  def auc(model, dataset)
    auc_by_model_and_dataset[model][dataset]
  end

  def dataset_quality(dataset)
    return nil  if models.empty? # it's impossible to calculate dataset weight without good models

    unless @dataset_qualities_cache && @dataset_qualities_cache[dataset] # caching
      @dataset_qualities_cache ||= {}
      @dataset_qualities_cache[dataset] = mean(aucs_for_dataset(dataset).values)
    end
    @dataset_qualities_cache[dataset]
  end

  # It's possible to calculate weighted AUC even for models which were excluded from @models
  # Weights of datasets though won't be recalculated and
  #   will consider only "good" models included in @models
  def weighted_auc(model, &dataset_quality_block)
    return nil  if datasets.empty?
    return nil  if models.empty? # In this case we can't estimate dataset weights
    dataset_quality_block = ->(dataset){ dataset_quality(dataset) } unless block_given?

    quality_norm_factor = datasets.map{|dataset|
      dataset_quality_block.call(dataset)
    }.inject(0.0, &:+)

    weighted_auc_total = datasets.map{|dataset|
      auc(model, dataset) * dataset_quality_block.call(dataset)
    }.inject(0.0, &:+)

    weighted_auc_total / quality_norm_factor
  end

  def best_model_among_collections(collections, banned_models: [])
    models.select{|model|
      collections.include?(model.collection_short_name)
    }.reject{|model|
      banned_models.include?(model)
    }.max_by{|model|
      weighted_auc(model)
    }
  end

  def without_bad_datasets(min_auc)
    return self.class.new(auc_by_model_and_dataset, models: [], datasets: [])  if models.empty?

    datasets_to_retain = datasets.select{|dataset|
      dataset_quality(dataset) >= min_auc
    }
    self.class.new(auc_by_model_and_dataset, models: models, datasets: datasets_to_retain)
  end

  def without_bad_models(min_auc)
    if datasets.empty?
      models_to_retain = []
    else
      models_to_retain = models.select{|model|
        weighted_auc(model) >= min_auc
      }
    end
    self.class.new(auc_by_model_and_dataset, models: models_to_retain, datasets: datasets)
  end

  def filter_datasets_by_species(species)
    AUCs.new(auc_by_model_and_dataset, models: models, datasets: datasets.select{|ds| ds.match(/_#{species}\./) })
  end

  def has_only_hocomoco_models?
    models.all?{|model| model.full_name.match(/~(DI)?HL~/) }
  end

  def has_chipseq_models?(uniprot, model_type)
    chipseq_dataset_code = (model_type == 'mono') ? 'CM' : 'CD'
    models.any?{|model| model.full_name.start_with?("#{uniprot}~#{chipseq_dataset_code}~") }
  end

  # Loads data for a single TF
  def self.from_folder(glob)
    auc_by_model_and_dataset = FileList[glob].map{|fn|
      model = Model.new_by_name(fn.pathmap('%n'))
      auc_by_dataset = File.readlines(fn).map{|line|
        dataset, auc, logauc = line.chomp.split("\t")
        [dataset, logauc.to_f]
      }.to_h
      [model, auc_by_dataset]
    }.to_h
    self.new(auc_by_model_and_dataset)
  end

  # Load hash {uniprot => AUCs} with iteratively filtered datasets and models
  def self.load_auc_infos_for_uniprot(min_weight_for_dataset: 0, min_auc_for_model: 0)
    semiuniprots = Dir.glob('auc/*').map{|fn| File.basename(fn).split('~')[0].split('_')[0] }.uniq

    result = semiuniprots.map{|uniprot|
      [uniprot, self.from_folder("auc/#{uniprot}_*.txt")]
    }.map{|uniprot, auc_infos|
      auc_infos_prev = nil
      while auc_infos != auc_infos_prev
        auc_infos_prev = auc_infos
        auc_infos = auc_infos.without_bad_datasets(min_weight_for_dataset).without_bad_models(min_auc_for_model)
      end
      [uniprot, auc_infos]
    }.to_h

    result.default_proc = ->(hsh,k) {
      self.new(Hash.new{|hsh2, k2| {}}, models: [], datasets: [])
    }
    result
  end

  def self.all_logaucs_in_folder(glob)
    initial_hsh = Hash.new{|h,k| h[k] = {} }
    FileList[glob].map{|fn|
      model = Model.new_by_name(fn.pathmap('%n'))
      auc_by_dataset = File.readlines(fn).map{|line|
        dataset, _auc, logauc = line.chomp.split("\t")
        [dataset, logauc.to_f]
      }.to_h
      [model, auc_by_dataset]
    } #.to_h
    .each_with_object(initial_hsh){|(model, auc_by_dataset), hsh|
      hsh[model].merge!(auc_by_dataset)
    }
  end

  def self.all_aucs_in_folder(glob)
    initial_hsh = Hash.new{|h,k| h[k] = {} }
    FileList[glob].map{|fn|
      model = Model.new_by_name(fn.pathmap('%n'))
      auc_by_dataset = File.readlines(fn).map{|line|
        dataset, auc, _logauc = line.chomp.split("\t")
        [dataset, auc.to_f]
      }.to_h
      [model, auc_by_dataset]
    } #.to_h
    .each_with_object(initial_hsh){|(model, auc_by_dataset), hsh|
      hsh[model].merge!(auc_by_dataset)
    }
  end

  def self.auc_infos_for_slice(all_aucs, slice_fn, model_type)
    case model_type.to_sym
    when :mono
      auc_infos_for_mono_slice(all_aucs, slice_fn)
    when :di
      auc_infos_for_di_slice(all_aucs, slice_fn)
    else
      raise "Undefined model type #{model_type}"
    end
  end

  def self.auc_infos_for_mono_slice(all_aucs, slice_fn)
    motif_sets = motifs_for_slice(slice_fn)
    models = []
    models += motif_sets[:chipseq_motifs].map{|motif_name|
      uniprot = motif_name.split('.').first
      Model.new("#{uniprot}~CM~#{motif_name}", :mono)
    }
    models += motif_sets[:hocomoco_motifs].map{|motif_name|
      uniprot = motif_name.split('.').first
      Model.new("#{uniprot}~HL~#{motif_name}", :mono)
    }
    models += motif_sets[:additional_motifs].map{|motif_name|
      uniprot = motif_name.split('.').first
      Model.new("#{uniprot}~AD~#{motif_name}", :mono)
    }

    result = models.select{|model| all_aucs.has_key?(model) }.map{|model|
      [model, all_aucs[model]]
    }.to_h
    self.new(result)
  end

  def self.auc_infos_for_di_slice(all_aucs, slice_fn)
    motif_sets = motifs_for_di_slice(slice_fn)
    models = []
    models += motif_sets[:chipseq_motifs].map{|motif_name|
      uniprot = motif_name.split('.').first
      Model.new("#{uniprot}~CD~#{motif_name}", :di)
    }
    models += motif_sets[:hocomoco_motifs].map{|motif_name|
      uniprot = motif_name.split('.').first
      Model.new("#{uniprot}~DIHL~#{motif_name}", :di)
    }
    models += motif_sets[:additional_motifs].map{|motif_name|
      uniprot = motif_name.split('.').first
      Model.new("#{uniprot}~DIAD~#{motif_name}", :di)
    }
    result = models.select{|model| all_aucs.has_key?(model) }.map{|model|
      [model, all_aucs[model]]
    }.to_h
    self.new(result)
  end

  def refined(min_weight_for_dataset: 0, min_auc_for_model: 0)
    auc_infos = self
    auc_infos_prev = nil
    while auc_infos != auc_infos_prev
      auc_infos_prev = auc_infos
      auc_infos = auc_infos.without_bad_datasets(min_weight_for_dataset).without_bad_models(min_auc_for_model)
    end
    auc_infos
  end
end
