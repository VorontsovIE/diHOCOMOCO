require 'models'
require 'best_models'

task :collection_coverage do
  auc_infos_for_uniprot = AUCs.load_auc_infos_for_uniprot(min_weight_for_dataset: 0.65,
                                                          min_auc_for_model: 0.65)

  secondary_models = File.readlines('secondary_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }
  banned_models = File.readlines('banned_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }
  best_models = collect_best_models(auc_infos_for_uniprot,
                                    secondary_models: secondary_models,
                                    banned_models: banned_models)

  models_by_collection = (Models::MonoCollections + Models::DiCollections).map{|collection|
    collection_models = Models.all_models.select{|model| model.collection_short_name == collection }
    [collection, collection_models]
  }.to_h.merge({'HOCOMOCO-v10' => best_models})

  ['HUMAN', 'MOUSE'].each do |species|
    ['mono', 'di'].each do |arity|
      count_tfs_by_collection = models_by_collection.map{|collection, models|
        [collection, models.select{|model| model.species == species && model.arity_type == arity }]
      }.map{|collection, models|
        [collection, models.map(&:uniprot).uniq.size]
      }.reject{|k,v| v.zero? }.sort_by{|k,v| k }.to_h
      
      puts '-----------------------------'
      puts "#{species}, #{arity}"
      count_tfs_by_collection.each{|collection, count|
        puts "#{Models::CollectionNames.fetch(collection){|k| k} }\t#{count}"
      }
    end
  end
end
