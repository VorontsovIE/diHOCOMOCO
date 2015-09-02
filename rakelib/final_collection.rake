require 'models'
require 'best_models'
require 'auc_info'
require 'auc_infos_filtering'

desc 'Collect final collection'
task :make_final_collection do
  filtering = AUCInfosFiltering.new(AUCInfo.load_all_infos)
  filtering.remove_bad_datasets_and_models!(0.6)

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
