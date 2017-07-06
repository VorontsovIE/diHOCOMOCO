require 'ape_precalculate_thresholds'

desc 'Precalculate thresholds'
task :precalculate_thresholds => [:precalculate_thresholds_mono, :precalculate_thresholds_di]

desc 'Precalculate thresholds for mononucleotide models'
task :precalculate_thresholds_mono

desc 'Precalculate thresholds for dinucleotide models'
task :precalculate_thresholds_di

task "precalculate_thresholds_mono" do
  SequenceDataset.each_dataset.map(&:uniprot).uniq.each do |uniprot|
    # Species defines background (it needn't match motif uniprot)
    ['HUMAN', 'MOUSE'].each do |species|
      pwm_folder = File.join('models/pwm/mono/all/', uniprot)
      next  unless Dir.exist?(pwm_folder)
      semiuniprot = uniprot.split('_').first
      background_fn = "control/local_backgrounds_averaged/#{semiuniprot}_#{species}.txt"
      next unless File.exist?(background_fn)

      output_folder = File.join('models/thresholds/mono/', "#{species}_background", uniprot)
      next  if Dir.exist?(output_folder) && FileList[File.join(pwm_folder, '*.pwm')].pathmap('%n').sort == FileList[File.join(output_folder, '*.thr')].pathmap('%n').sort
      FileUtils.rm_rf output_folder
      FileUtils.mkdir_p output_folder
      puts Ape.precalculate_thresholds_cmd pwm_folder,
                                          output_folder: output_folder,
                                          background: File.read(background_fn), # we always use dinucleotide background
                                          threshold_grid: ['1e-15', '1.0', '1.01', 'mul'],
                                          discretization: 1000,
                                          additional_options: ['--from-mono'], # same reasons: dinucleotide background
                                          mode: :di
    end
  end
end

=begin
task "precalculate_thresholds_di" do
  SequenceDataset.each_dataset do |control|
    pwm_folder = File.join('models/pwm_for_curation/di/all/', control.uniprot)
    output_folder = File.join('models/thresholds/di/', control.name)
    next  if Dir.exist?(output_folder)  &&  FileList[File.join(pwm_folder, '*.dpwm')].pathmap('%n').sort == FileList[File.join(output_folder, '*.thr')].pathmap('%n').sort
    rm_rf output_folder
    puts Ape.precalculate_thresholds_cmd pwm_folder,
                                         output_folder: output_folder,
                                         background: File.read(control.local_di_background_path),
                                         threshold_grid: ['1e-15', '1.0', '1.01', 'mul'],
                                         discretization: 1000,
                                         mode: :di
  end
end
=end
