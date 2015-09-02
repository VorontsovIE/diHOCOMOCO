require 'models'

def choose_best_collection(model_infos_by_collection, collections_of_interest)
  collections_of_interest.max_by{|collection|
    model_name, auc = model_infos_by_collection[collection]
    auc || 0
  }
end

def best_model_infos(model_infos_by_collection, quality_assessor, collections_of_interest)
  if (collections_of_interest & model_infos_by_collection.keys).empty?
    return { collection: nil, model: nil, auc: nil, quality: nil }
  end

  best_collection = choose_best_collection(model_infos_by_collection, collections_of_interest)
  best_model, best_auc = model_infos_by_collection[best_collection]
  {
    collection: best_collection,
    collection_fullname: Models::CollectionNames[best_collection],
    model: best_model,
    auc: best_auc,
    quality: quality_assessor.calculate_quality(best_model)
  }
end


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
