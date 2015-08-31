require 'models'
require 'sarus_results'

desc 'Calculate scores for each model on each control'
task :calculate_occurence_scores => [:calculate_occurence_scores_mono, :calculate_occurence_scores_di]

desc 'Calculate scores for each mononucleotide model on each control'
task :calculate_occurence_scores_mono

desc 'Calculate scores for each dinucleotide model on each control'
task :calculate_occurence_scores_di

SequenceDataset.each_dataset do |control|
  task "calculate_occurence_scores_mono:#{control.name}" do
    control.mono_models.each do |model|
      output_file = File.join('occurences/scores/mono/', control.uniprot, model.full_name, "#{control.name}.txt")
      next  if File.exist?(output_file)
      Sarus.run_besthits  control.filename,
                          model.path_to_pwm,
                          output_file: output_file,
                          mode: :mono
    end
  end
  task :calculate_occurence_scores_mono => "calculate_occurence_scores_mono:#{control.name}"

  task "calculate_occurence_scores_di:#{control.name}" do
    control.di_models.each do |model|
      output_file = File.join('occurences/scores/di/', control.uniprot, model.full_name, "#{control.name}.txt")
      next  if File.exist?(output_file)
      Sarus.run_besthits  control.filename,
                          model.path_to_pwm,
                          output_file: output_file,
                          mode: :di
    end
  end
  task :calculate_occurence_scores_di => "calculate_occurence_scores_di:#{control.name}"
end
