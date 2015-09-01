require 'ape_precalculate_thresholds'

desc 'Precalculate thresholds'
task :precalculate_thresholds => [:precalculate_thresholds_mono, :precalculate_thresholds_di]

desc 'Precalculate thresholds for mononucleotide models'
task :precalculate_thresholds_mono

desc 'Precalculate thresholds for dinucleotide models'
task :precalculate_thresholds_di

SequenceDataset.each_dataset do |control|
  task "precalculate_thresholds_mono:#{control.name}" do
    pwm_folder = File.join('models/pwm/mono/all/', control.uniprot)
    output_folder = File.join('models/thresholds/mono/all/', control.name)
    next  if Dir.exist?(output_folder) && FileList[File.join(pwm_folder, '*.pwm')].pathmap('%n').sort == FileList[File.join(output_folder, '*.thr')].pathmap('%n').sort
    Ape.run_precalculate_thresholds pwm_folder,
                                    output_folder: output_folder,
                                    background: File.read(control.local_di_background_path), # we always use dinucleotide background
                                    threshold_grid: ['1e-15', '1.0', '1.01', 'mul'],
                                    discretization: 1000,
                                    additional_options: ['--from-mono'], # same reasons: dinucleotide background
                                    mode: :di
  end
  task :precalculate_thresholds_mono => "precalculate_thresholds_mono:#{control.name}"


  task "precalculate_thresholds_di:#{control.name}" do
    pwm_folder = File.join('models/pwm/di/all/', control.uniprot)
    output_folder = File.join('models/thresholds/di/all/', control.name)
    next  if Dir.exist?(output_folder)  &&  FileList[File.join(pwm_folder, '*.dpwm')].pathmap('%n').sort == FileList[File.join(output_folder, '*.thr')].pathmap('%n').sort
    Ape.run_precalculate_thresholds pwm_folder,
                                    output_folder: output_folder,
                                    background: File.read(control.local_di_background_path),
                                    threshold_grid: ['1e-15', '1.0', '1.01', 'mul'],
                                    discretization: 1000,
                                    mode: :di
  end
  task :precalculate_thresholds_di => "precalculate_thresholds_di:#{control.name}"
end
