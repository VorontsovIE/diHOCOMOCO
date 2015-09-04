desc 'Remove models, failed curration'
task :remove_non_currated_models do
  # Datasets not in check_result.tsv are treated as good and passing curration
  # (thus corresponding models are not removed).
  #
  # There're such models and datasets from Papatsenko et al. (2015)
  # "Single-Cell Analyses of ESCs Reveal Alternative Pluripotent Cell States and
  #  Molecular Mechanisms that Control Self-Renewal."
  dataset_qualities = DatasetQuality.each_in_tsv('check_result.tsv').to_a
  models_failed_curration = dataset_qualities.reject(&:pass_quality_control?).flat_map(&:models)
  models_failed_curration.each do |model|
    $stderr.puts "Remove non-currated model #{model.full_name}"
    rm_f model.path_to_pcm
    rm_f model.path_to_pwm

    # Normally done above is enough but
    # if one has already calculated thresholds etc
    # we'd better remove them too

    SequenceDataset.each_for_uniprot(model.uniprot) do |control|
      if model.arity_type == 'mono'
        rm_f File.join('models/thresholds/mono/all/', control.name, "#{model.full_name}.thr")
      elsif model.arity_type == 'di'
        rm_f File.join('models/thresholds/di/all/', control.name, "#{model.full_name}.thr")
      else
        raise 'Unknown model arity'
      end
    end

    rm_rf  File.join('occurences/scores/mono/', model.uniprot, model.full_name)
    rm_rf  File.join('occurences/scores/di/', model.uniprot, model.full_name)
    rm_rf  File.join('occurences/pvalues/', model.uniprot, model.full_name)
    rm_rf  File.join('occurences/corrected_pvalues/', model.uniprot, model.full_name)
  end
end
