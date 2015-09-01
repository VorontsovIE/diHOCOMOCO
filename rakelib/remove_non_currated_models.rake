desc 'Remove models, failed curration'
task :remove_non_currated_models do
  dataset_qualities = DatasetQuality.each_in_xlsx('check_result.xlsx').to_a
  models_failed_curration = dataset_qualities.reject(&:pass_quality_control?).flat_map(&:models)
  models_failed_curration.each do |model|
    $stderr.puts "Remove non-currated model #{model.full_name}"
    rm model.path_to_pcm
    rm model.path_to_pwm

    SequenceDataset.each_for_uniprot(model.uniprot) do |control|
      rm File.join('models/thresholds/mono/all/', control.name, "#{model.full_name}.thr")
      rm File.join('models/thresholds/di/all/', control.name, "#{model.full_name}.thr")
    end

    rm_rf  File.join('occurences/scores/mono/', model.uniprot, model.full_name)
    rm_rf  File.join('occurences/scores/di/', model.uniprot, model.full_name)
    rm_rf  File.join('occurences/pvalues/', model.uniprot, model.full_name)
    rm_rf  File.join('occurences/corrected_pvalues/', model.uniprot, model.full_name)
  end
end
