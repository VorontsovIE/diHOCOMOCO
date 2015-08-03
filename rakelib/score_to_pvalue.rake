require 'sarus_results'
require 'find_pvalue_results'

desc 'Convert scores of motif occurences into pvalues'
task 'scores_to_pvalues' => ['scores_to_pvalues_mono', 'scores_to_pvalues_di']

SequenceDataset.each_file_by_glob('control/control/*.mfa') do |control|
  di_cmd = ['java', '-Xmx1G', '-cp', 'ape-2.0.1.jar', 'ru.autosome.ape.di.FindPvalue']
  task "scores_to_pvalues_di:#{control.name}" do
    Models.di_models.each do |model|
      scores_fn = File.join('occurences/scores/di/', control.uniprot, model.full_name, "#{control.name}.txt")
      scores = SarusResults.each_in_file(scores_fn).map(&:score)
      precalc_opts = ['--precalc', File.join('models/thresholds/di/all/', control.name)]
      
      background_fn = control_fn.pathmap('control/local_backgrounds/di/%n.txt')
      background_opt = ['--background', File.read(background_fn)]
      
      output_fn = File.join('occurences/pvalues/di/', control.uniprot, model.full_name, "#{control.name}.txt")
      sh *di_cmd, model.path_to_pwm, *scores.map(&:to_s) *pvalues, *precalc_opts, *background_opt, {out: output_fn}, {} 
    end
  end
  task 'scores_to_pvalues_di' => "scores_to_pvalues_di:#{control.name}"
end

SequenceDataset.each_file_by_glob('control/control/*.mfa') do |control|
  di_cmd = ['java', '-Xmx1G', '-cp', 'ape-2.0.1.jar', 'ru.autosome.ape.di.FindPvalue']
  task "scores_to_pvalues_mono:#{control.name}" do
    Models.mono_models.each do |model|
      scores_fn = File.join('occurences/scores/mono/', control.uniprot, model.full_name, "#{control.name}.txt")
      scores = SarusResults.each_in_file(scores_fn).map(&:score)
      precalc_opts = ['--precalc', File.join('models/thresholds/mono/all/', control.name)]

      background_fn = control_fn.pathmap('control/local_backgrounds/di/%n.txt')
      background_opt = ['--background', File.read(background_fn)]

      output_fn = File.join('occurences/pvalues/mono/', control.uniprot, model.full_name, "#{control.name}.txt")
      sh *di_cmd, model.path_to_pwm, *scores.map(&:to_s) *pvalues, *precalc_opts, *background_opt, '--from-mono', {out: output_fn}, {} 
    end
  end
  task 'scores_to_pvalues_mono' => "scores_to_pvalues_mono:#{control.name}"
end
