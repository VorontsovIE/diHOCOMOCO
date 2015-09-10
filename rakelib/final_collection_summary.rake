require 'models'
require 'best_models'
require 'quality_assessor'

def best_model_infos(auc_infos, quality_assessor, collections_of_interest)
  best_model = auc_infos.best_model_among_collections(collections_of_interest)
  if !best_model
    return { collection: nil, collection_fullname: nil, model: nil, model_name: nil, auc: nil, quality: nil }
  end

  {
    collection: best_model.collection_short_name,
    collection_fullname: Models::CollectionNames[best_model.collection_short_name],
    model: best_model,
    model_name: best_model.full_name,
    auc: auc_infos.weighted_auc(best_model),
    quality: quality_assessor.calculate_quality(best_model)
  }
end


desc 'Make final collection summary table (which model won etc)'
task :final_collection_summary do
  auc_infos_for_uniprot = AUCs.load_auc_infos_for_uniprot(min_weight_for_dataset: 0.65,
                                                          min_auc_for_model: 0.65)

  secondary_models = File.readlines('secondary_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }
  banned_models = File.readlines('banned_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }
  best_models = collect_best_models(auc_infos_for_uniprot,
                                    secondary_models: secondary_models,
                                    banned_models: banned_models)

  quality_assessor = QualityAssessor.new(auc_infos_for_uniprot, best_models: best_models, secondary_models: secondary_models)

  collections = (Models::MonoCollections + Models::DiCollections).sort

  headers = [
    'Uniprot', 'Species',
    'wAUC of best model', 'Best collection', 'Best model', 'Best model quality',
    'Number of datasets for TF',
    'Number of hocomoco models for TF',
    'Mono or dinucleotide win',
    'wAUC of best mononucleotide model', 'Best mononucleotide collection', 'Best mononucleotide model', 'Best mononucleotide model quality',
    'wAUC of best dinucleotide model', 'Best dinucleotide collection', 'Best dinucleotide model', 'Best dinucleotide model quality',
    *collections,
    'H10MO model', 'H10MO collection', 'H10MO wAUC', 'H10MO quality',
    'H10DI model', 'H10DI collection', 'H10DI wAUC', 'H10DI quality'
  ]

  results = auc_infos_for_uniprot.select{|uniprot, auc_infos|
    auc_infos.has_validation?
  }.map{|uniprot, auc_infos|
    auc_infos = auc_infos_for_uniprot[uniprot]
    best_infos_total = best_model_infos(auc_infos, quality_assessor, Models::MonoCollections + Models::DiCollections)
    best_infos_mono = best_model_infos(auc_infos, quality_assessor, Models::MonoCollections)
    best_infos_di = best_model_infos(auc_infos, quality_assessor, Models::DiCollections)

    num_hocomoco_models = Models.mono_models_by_uniprot(uniprot).count{|model|
      model.collection_short_name == 'HL'
    }

    aucs = collections.map{|collection|
      best_model_infos(auc_infos, quality_assessor, [collection])[:auc]
    }

    ############
    final_mono_models = best_models.select{|model|
      model.uniprot == uniprot
    }.select{|model|
      model.arity_type == 'mono'
    }.reject{|model|
      quality_assessor.calculate_quality(model) == 'S'
    }

    raise "WTF: #{final_mono_models.join(', ')}" if final_mono_models.size > 1
    final_mono_model = final_mono_models.first

    if final_mono_model
      final_mono_model_name = final_mono_model.full_name
      final_mono_model_collection = Models::CollectionNames[final_mono_model.collection_short_name]
      final_mono_model_auc = auc_infos.weighted_auc(final_mono_model)
      final_mono_model_quality = quality_assessor.calculate_quality(final_mono_model)
    end

    ############
    final_di_models = best_models.select{|model|
      model.uniprot == uniprot
    }.select{|model|
      model.arity_type == 'di'
    }.reject{|model|
      quality_assessor.calculate_quality(model) == 'S'
    }

    raise 'WTF: #{final_di_models.join(', ')}' if final_di_models.size > 1
    final_di_model = final_di_models.first

    if final_di_model
      final_di_model_name = final_di_model.full_name
      final_di_model_collection = Models::CollectionNames[final_di_model.collection_short_name]
      final_di_model_auc = auc_infos.weighted_auc(final_di_model)
      final_di_model_quality = quality_assessor.calculate_quality(final_di_model)
    end

    ############

    [ uniprot, uniprot[/_(?<species>HUMAN|MOUSE)$/, :species],
      *best_infos_total.values_at(:auc, :collection_fullname, :model_name, :quality),
      auc_infos.datasets.size,
      num_hocomoco_models,
      Models::MonoCollections.include?(best_infos_total[:collection]) ? 'Mono' : 'Di',
      *best_infos_mono.values_at(:auc, :collection_fullname, :model_name, :quality),
      *best_infos_di.values_at(:auc, :collection_fullname, :model_name, :quality),
      *aucs,
      final_mono_model_name, final_mono_model_collection, final_mono_model_auc, final_mono_model_quality,
      final_di_model_name, final_di_model_collection, final_di_model_auc, final_di_model_quality,
    ]
  }

  File.open('collection_perfomances.tsv', 'w') do |fw|
    fw.puts headers.join("\t")
    results.each{|row|
      fw.puts row.join("\t")
    }
  end
end
