require 'models'

def collection_checkers(collections)
  collections.map{|collection_checker|
    if collection_checker.respond_to?(:call)
      collection_checker
    else
      ->(model){ model.collection_short_name == collection_checker }
    end
  }
end

def checker_index(checkers, obj)
  checkers.index{|checker| checker.call(obj) }
end

def most_reliable_models(models, reliable_collections_ordered, banned: [])
  checkers = collection_checkers(reliable_collections_ordered)
  min_checker_index = models.reject{|model|
    banned.include?(model)
  }.map{|model|
    checker_index(checkers, model)
  }.compact.min
  return []  if !min_checker_index

  models.select{|model| checker_index(checkers, model) == min_checker_index }
end
