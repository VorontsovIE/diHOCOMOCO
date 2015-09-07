require 'models'
require 'best_models'
require 'auc_infos'
require 'html_table_output'

desc 'Collect final collection'
task :make_final_collection do
  # We add all secondary models explicitly at the last step
  human_secondary = File.readlines('human_secondary.txt').map(&:chomp).map{|name| Model.new_by_name(name) }
  mouse_secondary = File.readlines('mouse_secondary.txt').map(&:chomp).map{|name| Model.new_by_name(name) }
  secondary_models = human_secondary + mouse_secondary

  # For banned
  #   model without validation - we take next model by reliability
  #   validated hocomoco model - we take next by AUC motif instead
  #   validated non-hocomoco model - we take corresponding hocomoco models instead
  #     (sometimes they are as good as chipseq models but are much shorter
  #     while have AUC less by 0.001 or so)
  banned = File.readlines('banned.txt').map(&:chomp).map{|name| Model.new_by_name(name) }
  banned_hocomoco = banned.select{|model| model.collection_short_name == 'HL' }
  banned_not_hocomoco = banned.select{|model| model.collection_short_name != 'HL' }

  uniprots = Models.all_uniprots
  auc_infos_for_uniprot = AUCs.load_auc_infos_for_uniprot(min_weight_for_dataset: 0.65, min_auc_for_model: 0.65)

  ['HUMAN', 'MOUSE'].each{|species|
    ['mono', 'di'].each{|arity|
      mkdir_p "final_bundle/#{species}/#{arity}/pcm/" unless Dir.exist?("final_bundle/#{species}/#{arity}/pcm/")
      mkdir_p "final_bundle/#{species}/#{arity}/pwm/" unless Dir.exist?("final_bundle/#{species}/#{arity}/pwm/")
      mkdir_p "final_bundle/#{species}/#{arity}/logo/" unless Dir.exist?("final_bundle/#{species}/#{arity}/logo/")
    }
  }

  best_models_mono = []
  best_models_di = []

  uniprots.each do |uniprot|
    auc_infos = auc_infos_for_uniprot[uniprot]
    best_model_mono = nil
    best_model_di = nil

    # If we have validated model, use it.
    # If we don't have, we take model based on our believes in collection qualities

    if auc_infos && !auc_infos.empty?
      # Model have been validated (but not necessary that both mono- and di- were)

      ### Mononucleotide models
      best_model_mono = auc_infos.best_model_among_collections(
        Models::MonoCollectionsForFinalBundle,
        banned: banned_hocomoco
      )
      if best_model_mono # Mono-model was validated
        if banned_not_hocomoco.include?(best_model_mono) # Validated model was a banned non-hocomoco model
          hocomoco_models = Models.mono_models_by_uniprot(uniprot).select{|model|
            model.collection_short_name == 'HL'
          }
          # So we add hocomoco models instead
          best_models_mono += hocomoco_models # both primary and secondary if they exist
          # If we added banned models here, we'll remove them manually at the last step
        else
          best_models_mono << best_model_mono
        end
      else # Mono-model wasn't validated
        best_models_mono += most_reliable_models(
                              Models.mono_models_by_uniprot(uniprot),
                              Models::MonoCollectionsReliability,
                              banned: banned
                            )
      end

      ### Dinucleotide models
      best_model_di = auc_infos.best_model_among_collections(
        Models::DiCollectionsForFinalBundle
      )
      if best_model_di # Dinucleotide model was validated
        # Only take di-model when it beats mono model or if mono model doesn't exist
        if best_model_mono
          if auc_infos.weighted_model_aucs[best_model_di] > auc_infos.weighted_model_aucs[best_model_mono]
            best_models_di << best_model_di
          end
        else
          best_models_di << best_model_di
        end
      else # Non-validate dinucleotide model
        best_models_di += most_reliable_models(
                            Models.di_models_by_uniprot(uniprot),
                            Models::DiCollectionsReliability,
                            banned: banned
                          )
      end
    else
      best_models_mono += most_reliable_models(
                            Models.mono_models_by_uniprot(uniprot),
                            Models::MonoCollectionsReliability,
                            banned: banned
                          )
      best_models_di += most_reliable_models(
                          Models.di_models_by_uniprot(uniprot),
                          Models::DiCollectionsReliability,
                          banned: banned
                        )
    end
  end

  # Also we take all secondary models
  best_models_mono += secondary_models.select{|model| model.mono_or_di_mode.arity_type == 'mono' }
  best_models_di += secondary_models.select{|model| model.mono_or_di_mode.arity_type == 'di' }

  best_models_mono = best_models_mono.compact.uniq
  best_models_di = best_models_di.compact.uniq

  # We don't expect banned models here but they may appear
  #   when we append hocomoco models instead of a banned one.
  # So it's more safe to remove them manually
  best_models_mono.reject!{|model| banned.include?(model) }
  best_models_di.reject!{|model| banned.include?(model) }

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
end
