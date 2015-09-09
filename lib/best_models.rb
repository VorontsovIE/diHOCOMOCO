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

def most_reliable_models(models, reliable_collections_ordered, banned_models: [])
  checkers = collection_checkers(reliable_collections_ordered)
  min_checker_index = models.reject{|model|
    banned_models.include?(model)
  }.map{|model|
    checker_index(checkers, model)
  }.compact.min
  return []  if !min_checker_index

  models.select{|model| checker_index(checkers, model) == min_checker_index }
end

# Collect all models for final collection
def collect_best_models(auc_infos_for_uniprot, secondary_models:, banned_models:)
  # We add all secondary models explicitly at the last step

  # For banned
  #   model without validation - we take next model by reliability
  #   validated hocomoco model - we take next by AUC motif instead
  #   validated non-hocomoco model - we take corresponding hocomoco models instead
  #     (sometimes they are as good as chipseq models but are much shorter
  #     while have AUC less by 0.001 or so)
  banned_hocomoco = banned_models.select{|model| model.collection_short_name == 'HL' }
  banned_not_hocomoco = banned_models.select{|model| model.collection_short_name != 'HL' }

  uniprots = Models.all_uniprots

  best_models_mono = []
  best_models_di = []

  uniprots.each do |uniprot|
    auc_infos = auc_infos_for_uniprot[uniprot]
    best_model_mono = nil
    best_model_di = nil

    # If we have validated model, use it.
    # If we don't have, we take model based on our believes in collection qualities

    if auc_infos.has_validation?
      # Model have been validated (but not necessary that both mono- and di- were)

      ### Mononucleotide models
      best_model_mono = auc_infos.best_model_among_collections(
        Models::MonoCollectionsForFinalBundle,
        banned_models: banned_hocomoco
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
                              banned_models: banned_models
                            )
      end

      ### Dinucleotide models
      best_model_di = auc_infos.best_model_among_collections(
        Models::DiCollectionsForFinalBundle
      )
      if best_model_di # Dinucleotide model was validated
        # Only take di-model when it beats mono model or if mono model doesn't exist
        if best_model_mono
          if auc_infos.weighted_auc(best_model_di) > auc_infos.weighted_auc(best_model_mono)
            best_models_di << best_model_di
          end
        else
          best_models_di << best_model_di
        end
      else # Non-validate dinucleotide model
        best_models_di += most_reliable_models(
                            Models.di_models_by_uniprot(uniprot),
                            Models::DiCollectionsReliability,
                            banned_models: banned_models
                          )
      end
    else
      best_models_mono += most_reliable_models(
                            Models.mono_models_by_uniprot(uniprot),
                            Models::MonoCollectionsReliability,
                            banned_models: banned_models
                          )
      best_models_di += most_reliable_models(
                          Models.di_models_by_uniprot(uniprot),
                          Models::DiCollectionsReliability,
                          banned_models: banned_models
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
  best_models_mono.reject!{|model| banned_models.include?(model) }
  best_models_di.reject!{|model| banned_models.include?(model) }

  best_models_mono + best_models_di
end


def same_by?(models, &block)
  characteristics = models.map(&block)
  characteristics.all?{|ch| ch == characteristics.first }
end

# Calculate a characteristic for model. Make sure that all models have the same value or raise.
def take_and_check_consistency(models, &block)
  raise 'Can\'t take characteristic for empty model list'  if models.empty?
  if same_by?(models, &block)
    block.call(models.first)
  else
    raise 'Inconsistent characteristics for joint models #{models.inspect}'
  end
end
