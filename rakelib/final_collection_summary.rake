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

  quality_assessor = QualityAssessor.new(auc_infos_for_uniprot)

  collections = (Models::MonoCollections + Models::DiCollections).sort

  headers = [
    'Uniprot', 'Species',
    'Max AUC of best model', 'Best collection', 'Best model', 'Best model quality',
    'Number of hocomoco models for TF',
    'Mono or dinucleotide win',
    'Max AUC of best mononucleotide model', 'Best mononucleotide collection', 'Best mononucleotide model', 'Best mononucleotide model quality',
    'Max AUC of best dinucleotide model', 'Best dinucleotide collection', 'Best dinucleotide model', 'Best dinucleotide model quality',
    *collections
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

    [ uniprot, uniprot[/_(?<species>HUMAN|MOUSE)$/, :species],
      *best_infos_total.values_at(:auc, :collection_fullname, :model_name, :quality),
      num_hocomoco_models,
      Models::MonoCollections.include?(best_infos_total[:collection]) ? 'Mono' : 'Di',
      *best_infos_mono.values_at(:auc, :collection_fullname, :model_name, :quality),
      *best_infos_di.values_at(:auc, :collection_fullname, :model_name, :quality),
      *aucs ]
  }

  File.open('collection_perfomances.tsv', 'w') do |fw|
    fw.puts headers.join("\t")
    results.each{|row|
      fw.puts row.join("\t")
    }
  end
end
