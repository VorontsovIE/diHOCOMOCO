desc 'Calculate scores for each model on each control'
task :calculate_occurence_scores => [:calculate_occurence_scores_mono, :calculate_occurence_scores_di]


desc 'Calculate scores for each model on each mononucleotide control'
task :calculate_occurence_scores_mono

FileList['control/control/*.mfa'].sort.each do |control_fn|
  cmd = ['java', '-cp', 'sarus.jar', 'ru.autosome.SARUS']
  control_name = control_fn.pathmap('%n')
  uniprot = control_name[/^.+_(HUMAN|MOUSE)/]
  output_dir = File.join('occurences/scores/mono/', uniprot)

  directory output_dir do
    mkdir_p output_dir

    FileList["models/pwm/mono/all/#{uniprot}/*.pwm"].each do |pwm_fn|
      pwm_name = pwm_fn.pathmap('%n')
      mkdir_p File.join(output_dir, pwm_name)
      output_fn = File.join(output_dir, pwm_name, "#{control_name}.txt")
      sh *cmd, control_fn, pwm_fn, 'besthit', 'suppress', {out: output_fn}, {}
    end
  end
  task :calculate_occurence_scores_mono => output_dir
end


desc 'Calculate scores for each model on each dinucleotide control'
task :calculate_occurence_scores_di

FileList['control/control/*.mfa'].sort.each do |control_fn|
  cmd = ['java', '-cp', 'sarus.jar', 'ru.autosome.di.SARUS']
  control_name = control_fn.pathmap('%n')
  uniprot = control_name[/^.+_(HUMAN|MOUSE)/]
  output_dir = File.join('occurences/scores/di/', uniprot)

  directory output_dir do
    mkdir_p output_dir

    FileList["models/pwm/di/all/#{uniprot}/*.dpwm"].each do |pwm_fn|
      pwm_name = pwm_fn.pathmap('%n')
      mkdir_p File.join(output_dir, pwm_name)
      output_fn = File.join(output_dir, pwm_name, "#{control_name}.txt")
      sh *cmd, control_fn, pwm_fn, 'besthit', 'suppress', {out: output_fn}, {}
    end
  end
  task :calculate_occurence_scores_di => output_dir
end

