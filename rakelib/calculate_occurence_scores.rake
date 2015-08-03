require 'sarus_results'

desc 'Calculate scores for each model on each control'
task :calculate_occurence_scores => [:calculate_occurence_scores_mono, :calculate_occurence_scores_di]

desc 'Calculate scores for each model on each mononucleotide control'
task :calculate_occurence_scores_mono

desc 'Calculate scores for each model on each dinucleotide control'
task :calculate_occurence_scores_di

SequenceDataset.each_file_by_glob('control/control/*.mfa') do |control|
  task "calculate_occurence_scores_mono:#{control.name}" do
    Models.mono_models.each do |model|
      Sarus.run_besthits  control.filename,
                          model.path_to_pwm,
                          output_file: File.join('occurences/scores/mono/', control.uniprot, model.full_name, "#{control.name}.txt"),
                          mode: :mono
    end
  end
  task :calculate_occurence_scores_mono => "calculate_occurence_scores_mono:#{control.name}"

  task "calculate_occurence_scores_di:#{control.name}" do
    Models.di_models.each do |model|
      Sarus.run_besthits  control.filename,
                          model.path_to_pwm,
                          output_file: File.join('occurences/scores/di/', control.uniprot, model.full_name, "#{control.name}.txt"),
                          mode: :di
    end
  end
  task :calculate_occurence_scores_di => "calculate_occurence_scores_di:#{control.name}"
end
