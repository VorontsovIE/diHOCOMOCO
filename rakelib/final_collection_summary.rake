require 'models'
require 'best_models'
require 'auc_info'
require 'auc_infos_filtering'
require 'quality_assessor'

desc 'Make final collection summary table (which model won etc)'
task :final_collection_summary do
  filtering = AUCInfosFiltering.new(AUCInfo.load_all_infos)
  filtering.remove_bad_datasets_and_models!(0.6)
  quality_assessor = QualityAssessor.new(filtering)
  # not_filtered_aucs = AUCInfo.load_all_infos # Don't use it in AUCInfosFiltering because it will modify this variable
  collection_perfomances = filtering.model_perfomances_collections_grouped

  collections = (Models::MonoCollections + Models::DiCollections).sort

  headers = [
    'Uniprot', 'Species',
    'Median AUC of best model', 'Best collection', 'Best model', 'Best model quality',
    'Number of hocomoco models for TF',
    'Mono or dinucleotide win',
    'Median AUC of best mononucleotide model', 'Best mononucleotide collection', 'Best mononucleotide model', 'Best mononucleotide model quality',
    'Median AUC of best dinucleotide model', 'Best dinucleotide collection', 'Best dinucleotide model', 'Best dinucleotide model quality',
    *collections
  ]

  results = []
  collection_perfomances.each do |uniprot, model_infos_by_collection|
    best_infos_total = best_model_infos(model_infos_by_collection, quality_assessor, collections)
    best_infos_mono = best_model_infos(model_infos_by_collection, quality_assessor, Models::MonoCollections)
    best_infos_di = best_model_infos(model_infos_by_collection, quality_assessor, Models::DiCollections)

    num_hocomoco_models = Models.mono_models_by_uniprot(uniprot).count{|model| model.collection_short_name == 'HL' }

    aucs = collections.map{|collection|
      model_name, auc = model_infos_by_collection[collection]
      auc
    }

    results << [
      uniprot, uniprot[/_(?<species>HUMAN|MOUSE)$/,:species],
      *best_infos_total.values_at(:auc, :collection, :model, :quality),
      num_hocomoco_models,
      Models::MonoCollections.include?(best_infos_total[:collection]) ? 'Mono' : 'Di',
      *best_infos_mono.values_at(:auc, :collection, :model, :quality),
      *best_infos_di.values_at(:auc, :collection, :model, :quality),
      *aucs
    ]
  end

  hocomoco_qualities = File.readlines('hocomoco_qualities.tsv').map{|line| line.chomp.split("\t") }.to_h

  # We don't add any dinucleotide models w/o validation, only mononucleotide
  best_models_wo_validation = best_models_wo_benchmark(collection_perfomances, ['HL', 'SMI', 'CM'], :mono)
  best_models_wo_validation.each do |best_model|
    model_collection = best_model.collection_short_name
    num_hocomoco_models = Models.mono_models_by_uniprot(best_model.uniprot).count{|model| model.collection_short_name == 'HL' }
    if model_collection == 'HL'
      quality = hocomoco_qualities[best_model.model_name]
    elsif model_collection == 'SMI'
      quality = 'E'
    else
      raise "Model `#{best_model}` can\'t be selected without benchmark. Only HL and SMI models can, CM can't be here too"
    end
    results << [
      best_model.uniprot, best_model.species,
      nil, model_collection, best_model.full_name, quality,
      num_hocomoco_models,
      'Mono',
      nil, model_collection, best_model.full_name, quality,
      nil, nil, nil, nil,
      *[nil] * collections.size
    ]
  end

  File.open('collection_perfomances.tsv', 'w') do |fw|
    fw.puts headers.join("\t")
    results.each{|row|
      fw.puts row.join("\t")
    }
  end
end
