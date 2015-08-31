require 'models'
require 'sarus_results'
require 'ape_find_pvalue'
require 'find_pvalue_results'

desc 'Convert scores of motif occurences into pvalues'
task :scores_to_pvalues => [:scores_to_pvalues_mono, :scores_to_pvalues_di]

desc 'Convert scores of motif occurences into pvalues for mononucleotide models'
task :scores_to_pvalues_mono

desc 'Convert scores of motif occurences into pvalues for dinucleotide models'
task :scores_to_pvalues_di

SequenceDataset.each_dataset do |control|
  task "scores_to_pvalues_mono:#{control.name}" do
    control.mono_models.each do |model|
      scores_fn = File.join('occurences/scores/mono/', control.uniprot, model.full_name, "#{control.name}.txt")
      scores = SarusResults.each_in_file(scores_fn).map(&:score)
      output_file = File.join('occurences/pvalues/', control.uniprot, model.full_name, "#{control.name}.txt")
      next  if File.exist?(output_file)
      Ape.run_find_pvalue   model.path_to_pwm,
                            scores,
                            output_file: output_file,
                            background: File.read(control.local_di_background_path),
                            discretization: 1000,
                            additional_options: ['--precalc', File.join('models/thresholds/mono/all/', control.name)] + ['--from-mono'],
                            mode: :di
    end
  end
  task :scores_to_pvalues_mono => "scores_to_pvalues_mono:#{control.name}"


  task "scores_to_pvalues_di:#{control.name}" do
    control.di_models.each do |model|
      scores_fn = File.join('occurences/scores/di/', control.uniprot, model.full_name, "#{control.name}.txt")
      scores = SarusResults.each_in_file(scores_fn).map(&:score)
      output_file = File.join('occurences/pvalues/', control.uniprot, model.full_name, "#{control.name}.txt")
      next  if File.exist?(output_file)
      Ape.run_find_pvalue   model.path_to_pwm,
                            scores,
                            output_file: output_file,
                            background: File.read(control.local_di_background_path),
                            discretization: 1000,
                            additional_options: ['--precalc', File.join('models/thresholds/di/all/', control.name)],
                            mode: :di
    end
  end
  task :scores_to_pvalues_di => "scores_to_pvalues_di:#{control.name}"
end
