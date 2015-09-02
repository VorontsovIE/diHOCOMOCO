require 'dataset_qualities'
require 'models'
require 'sequence_dataset'
require 'auc_info'
require 'auc_infos_filtering'
require 'quality_assessor'

# validated models are taken from collections of interest
def uniprots_with_validated_models(collection_perfomances, collections_of_interest)
  Models.all_uniprots.select{|uniprot|
    collection_perfomances.has_key?(uniprot)
  }.reject{|uniprot|
    (collection_perfomances[uniprot].keys & collections_of_interest).empty?
  }
end

# validated models are taken from collections of interest
def uniprots_without_validated_models(collection_perfomances, collections_of_interest)
  Models.all_uniprots - uniprots_with_validated_models(collection_perfomances, collections_of_interest)
end

def best_models_with_benchmark(collection_perfomances, collections_of_interest, arity_type)
  raise "Unknown arity type `#{arity_type}`. Should be :mono or :di"  unless [:mono, :di].include?(arity_type)
  uniprots_with_validated_models(collection_perfomances, collections_of_interest).map{|uniprot|
    collection_perfomances[uniprot].select{|collection, (model, auc)|
      collections_of_interest.include?(collection)
    }.map{|collection, (model, auc)|
      [model, auc]
    }.max_by{|model, auc|
      auc
    }.first
  }.map{|model_fullname|
    Model.new(model_fullname, arity_type)
  }
end

# ordered_collections_of_interest looks like ['HL', 'SMI', 'CM']
def best_models_wo_benchmark(collection_perfomances, ordered_collections_of_interest, arity_type)
  raise "Unknown arity type `#{arity_type}`. Should be :mono or :di"  unless [:mono, :di].include?(arity_type)

  uniprots_without_validated_models(collection_perfomances, ordered_collections_of_interest).map{|uniprot|
    models = (arity_type == :mono) ? Models.mono_models_by_uniprot(uniprot) : Models.di_models_by_uniprot(uniprot)
    models.select{|model|
      ordered_collections_of_interest.include?(model.collection_short_name)
    }
  }.reject(&:empty?).map{|models_for_uniprot|
    best_collection = models_for_uniprot.map(&:collection_short_name).min_by{|collection|
      ordered_collections_of_interest.index(collection)
    }
    models_for_uniprot.select{|model|
      model.collection_short_name == best_collection
    }
  }.flatten # Here for a single TF we can have several models w/o benchmark
end

desc 'Collect final collection'
task :make_final_collection do
  dataset_quality_by_name = DatasetQuality.each_in_xlsx('check_result.xlsx').map{|infos|
    [infos.dataset_name, infos]
  }.to_h

  filtering = AUCInfosFiltering.new(AUCInfo.load_all_infos)
  filtering.remove_bad_datasets_and_models!(dataset_quality_by_name, 0.6)

  collection_perfomances = filtering.model_perfomances_collections_grouped

  ['HUMAN', 'MOUSE'].each{|species|
    mkdir_p "final_bundle/#{species}/mono/pcm/"
    mkdir_p "final_bundle/#{species}/mono/pwm/"
  }
  best_models_mono = []
  best_models_mono += best_models_with_benchmark(collection_perfomances, Models::CollectionsForFinalBundle & Models::MonoCollections, :mono)
  best_models_mono += best_models_wo_benchmark(collection_perfomances, ['HL', 'SMI', 'CM'], :mono)
  best_models_mono.each{|model|
    cp model.path_to_pcm, "final_bundle/#{model.species}/mono/pcm/"
    cp model.path_to_pwm, "final_bundle/#{model.species}/mono/pwm/"
  }

  ['HUMAN', 'MOUSE'].each{|species|
    mkdir_p "final_bundle/#{species}/di/pcm/"
    mkdir_p "final_bundle/#{species}/di/pwm/"
  }
  best_models_di = []
  best_models_di += best_models_with_benchmark(collection_perfomances, Models::CollectionsForFinalBundle & Models::DiCollections, :di)
  # Don't add any dinucleotide models w/o benchmark
  # best_models_di += best_models_wo_benchmark(collection_perfomances, ['CD'], :di)
  best_models_di.each{|model|
    cp model.path_to_pcm, "final_bundle/#{model.species}/di/pcm/"
    cp model.path_to_pwm, "final_bundle/#{model.species}/di/pwm/"
  }
end
