require 'median'
require 'models'

class AUCs
  attr_reader :auc_by_model_and_dataset
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

  def models
    auc_by_model_and_dataset.keys
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

  # def best_model_among_collections(collections)
  #    best_model_mono, best_auc_mono = model_aucs.select{|model, auc|
  #     (Models::CollectionsForFinalBundle & Models::MonoCollections).include?(model.collection_short_name)
  #   }.max_by{|model, auc| auc }
  # end

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
    [uniprot, auc_infos.without_bad_datasets(0.65)]
  }.to_h;

  model_aucs_for_uniprot = auc_infos_for_uniprot.map{|uniprot, auc_infos|
    [uniprot, auc_infos.weighted_model_aucs]
  }.map{|uniprot, model_aucs|
    good_model_aucs = model_aucs.select{|model, auc| auc >= 0.65 }
    [uniprot, good_model_aucs]
  }.reject{|uniprot, model_aucs|
    model_aucs.empty?
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
  model_aucs_for_uniprot.each{|uniprot, model_aucs|
    best_model_mono, best_auc_mono = model_aucs.select{|model, auc|
      (Models::CollectionsForFinalBundle & Models::MonoCollections).include?(model.collection_short_name)
    }.max_by{|model, auc| auc }

    best_model_di, best_auc_di = model_aucs.select{|model, auc|
      (Models::CollectionsForFinalBundle & Models::DiCollections).include?(model.collection_short_name)
    }.max_by{|model, auc| auc }

    next  unless best_model_mono || best_model_di

    num_models += 1
    
    print '<tr class="mono">'
    print "<td rowspan=2>#{uniprot}</td>"
    if best_model_mono
      num_mono_models += 1
      print "<td>#{best_auc_mono}</td>"
      print "<td>#{best_model_mono.full_name}</td>"
      print "<td><img src='#{best_model_mono.path_to_logo}'/></td>"
    end
    print '</tr>'
    print '<tr class="di">'
    if best_model_di
      num_di_models += 1
      print "<td>#{best_auc_di}</td>"
      print "<td>#{best_model_di.full_name}</td>"
      print "<td><img src='#{best_model_di.path_to_logo}'/></td>"
    end
    puts "</tr>"
  }

  puts '</table>'
  puts '</body></html>'

  $stderr.puts({models: num_models, mono_models: num_mono_models, di_models: num_di_models})
end
