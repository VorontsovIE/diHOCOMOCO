# Наверное надо не по всем контролям, а только по контролям мотива
desc 'Precalculate thresholds'
task :precalculate_thresholds => [:precalculate_thresholds_mono, :precalculate_thresholds_di]


desc 'Precalculate thresholds for mononucleotide models'
task :precalculate_thresholds_mono

FileList['control/control/*'].each do |control_fn|
  control_name = control_fn.pathmap('%n')
  output_dir = File.join('models/thresholds/mono/all/', control_name)
  directory output_dir do
    mkdir_p output_dir
    uniprot = control_name[/^.+_(HUMAN|MOUSE)/]
    motif_dir = File.join('models/pwm/mono/all/', uniprot)

    background_fn = control_fn.pathmap('control/local_backgrounds/mono/%n.txt')
    background = File.read(background_fn).split.map(&:to_f)
    background_opt = ['--background', background.join(',')]

    script_cmd = ['java', '-Xmx1G', '-cp', 'ape-2.0.1.jar', 'ru.autosome.ape.PrecalculateThresholds']
    threshold_grid = ['--pvalues', ['1e-15', '1.0', '1.05', 'mul'].join(',')]
    sh *script_cmd, motif_dir, output_dir, *threshold_grid, '--silent', '--discretization', '1000', *background_opt
  end
  task :precalculate_thresholds_mono => output_dir
end


desc 'Precalculate thresholds for dinucleotide models'
task :precalculate_thresholds_di

FileList['control/control/*'].each do |control_fn|
  control_name = control_fn.pathmap('%n')
  output_dir = File.join('models/thresholds/di/all/', control_name)
  directory output_dir do
    mkdir_p output_dir
    uniprot = control_name[/^.+_(HUMAN|MOUSE)/]
    motif_dir = File.join('models/pwm/di/all/', uniprot)

    background_fn = control_fn.pathmap('control/local_backgrounds/di/%n.txt')
    background = File.read(background_fn).split.map(&:to_f)
    background_opt = ['--background', background.join(',')]

    script_cmd = ['java', '-Xmx1G', '-cp', 'ape-2.0.1.jar', 'ru.autosome.ape.di.PrecalculateThresholds']
    threshold_grid = ['--pvalues', ['1e-15', '1.0', '1.05', 'mul'].join(',')]
    sh *script_cmd, motif_dir, output_dir, *threshold_grid, '--silent', '--discretization', '1000', *background_opt
  end
  task :precalculate_thresholds_di => output_dir
end
