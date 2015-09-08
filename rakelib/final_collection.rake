require 'models'
require 'best_models'
require 'auc_infos'
require 'html_table_output'

desc 'Collect final collection'
task :make_final_collection do
  auc_infos_for_uniprot = AUCs.load_auc_infos_for_uniprot(min_weight_for_dataset: 0.65,
                                                          min_auc_for_model: 0.65)

  secondary_models = File.readlines('secondary_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }
  banned_models = File.readlines('banned_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }
  best_models_mono, best_models_di = collect_best_models( auc_infos_for_uniprot,
                                                          secondary_models: secondary_models,
                                                          banned_models: banned_models )
  best_models = best_models_mono + best_models_di

  quality_assessor = QualityAssessor.new(auc_infos_for_uniprot, secondary_models)

  File.open('final_collection.html', 'w') do |fw|
    print_html_table_for_grouped_models(
      auc_infos_for_uniprot,
      best_models.group_by(&:uniprot),
      quality_assessor,
      stream: fw
    )
  end

  ['HUMAN', 'MOUSE'].each{|species|
    ['mono', 'di'].each{|arity|
      mkdir_p "final_bundle/#{species}/#{arity}/pcm/" unless Dir.exist?("final_bundle/#{species}/#{arity}/pcm/")
      mkdir_p "final_bundle/#{species}/#{arity}/pwm/" unless Dir.exist?("final_bundle/#{species}/#{arity}/pwm/")
      mkdir_p "final_bundle/#{species}/#{arity}/logo/" unless Dir.exist?("final_bundle/#{species}/#{arity}/logo/")
    }
  }

  best_models_mono.each{|model|
    cp model.path_to_pcm, "final_bundle/#{model.species}/mono/pcm/"
    cp model.path_to_pwm, "final_bundle/#{model.species}/mono/pwm/"
    cp model.path_to_logo, "final_bundle/#{model.species}/mono/logo/"
  }

  best_models_di.each{|model|
    cp model.path_to_pcm, "final_bundle/#{model.species}/di/pcm/"
    cp model.path_to_pwm, "final_bundle/#{model.species}/di/pwm/"
    cp model.path_to_logo, "final_bundle/#{model.species}/di/logo/"
  }
end
