require 'median'
require 'models'

# AUCs represents AUC values for a single TF (several models, several datasets; AUC for each model over each dataset)
class AUCs
  attr_reader :auc_by_model_and_dataset

  # {model => {dataset => auc}}
  # model should be an instance of class Model (see #best_model_among_collections)
  def initialize(auc_by_model_and_dataset)
    @auc_by_model_and_dataset = auc_by_model_and_dataset
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
      }.to_h
    end
  end

  def best_model_among_collections(collections)
     best_model, best_auc = weighted_model_aucs.select{|model, auc|
      collections.include?(model.collection_short_name)
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
end


task :model_logos_table do
  uniprots = FileList['occurences/auc/*'].pathmap('%n');
  auc_infos_for_uniprot = uniprots.map{|uniprot|
    [uniprot, AUCs.from_folder("occurences/auc/#{uniprot}/*.txt")]
  }.map{|uniprot, auc_infos|
    auc_infos_prev = nil
    while auc_infos != auc_infos_prev
      auc_infos_prev = auc_infos
      auc_infos = auc_infos.without_bad_datasets(0.65).without_bad_models(0.65)
    end
    [uniprot, auc_infos]
  }.reject{|uniprot,auc_infos|
    auc_infos.empty?
  }.to_h;

  puts '<html><head><style>'
  puts 'table, tr, td{ border: 1px solid black; border-collapse: collapse; }'
  puts 'tr.mono, tr.mono td { border-top: 3px solid black; }'
  puts 'tr.di, tr.di td{ border-bottom: 3px solid black; }'
  puts 'tr.mono td:first-child{ border-bottom: 3px solid black; border-left: 3px solid;}'
  puts 'tr td:last-child{ border-right: 3px solid; }'
  puts '</style></head><body>'
  puts '<table>'

  num_models = num_mono_models = num_di_models = 0
  auc_infos_for_uniprot.each{|uniprot, auc_infos|
    best_model_mono = auc_infos.best_model_among_collections(Models::CollectionsForFinalBundle & Models::MonoCollections)
    best_model_di = auc_infos.best_model_among_collections(Models::CollectionsForFinalBundle & Models::DiCollections)

    next  unless best_model_mono || best_model_di

    num_models += 1

    print '<tr class="mono">'
    print "<td rowspan=2>#{uniprot}</td>"
    if best_model_mono
      num_mono_models += 1
      best_auc_mono = auc_infos.weighted_model_aucs[best_model_mono]
      print "<td>#{best_auc_mono.round(3)}</td>"
      print "<td>#{best_model_mono.full_name}</td>"
      print "<td><img src='#{best_model_mono.path_to_logo}'/></td>"
    end
    print '</tr>'
    print '<tr class="di">'
    if best_model_di
      num_di_models += 1
      best_auc_di = auc_infos.weighted_model_aucs[best_model_di]
      print "<td>#{best_auc_di.round(3)}</td>"
      print "<td>#{best_model_di.full_name}</td>"
      print "<td><img src='#{best_model_di.path_to_logo}'/></td>"
    end
    puts "</tr>"
  }

  puts '</table>'
  puts '</body></html>'

  $stderr.puts({models: num_models, mono_models: num_mono_models, di_models: num_di_models})
end
