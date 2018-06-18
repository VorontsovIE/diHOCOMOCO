require 'ape_precalculate_thresholds'

desc 'Precalculate thresholds'
task :precalculate_thresholds => [:precalculate_thresholds_mono, :precalculate_thresholds_di]

desc 'Precalculate thresholds for mononucleotide models'
task :precalculate_thresholds_mono

desc 'Precalculate thresholds for dinucleotide models'
task :precalculate_thresholds_di

['mono', 'di'].each do |model_type|
  task "precalculate_thresholds_#{model_type}" do
    additional_options = ['--single-motif']
    additional_options << '--from-mono'  if model_type == 'mono'  # we always use dinucleotide background
    SequenceDataset.each_dataset.map(&:uniprot).uniq.each do |uniprot|
      # Species defines background (it needn't match motif uniprot)
      ['HUMAN', 'MOUSE'].each do |species|
        pwm_folder = File.join("models/pwm/#{model_type}/all/", uniprot)
        next  unless Dir.exist?(pwm_folder)
        semiuniprot = uniprot.split('_').first
        background_fn = "control/local_backgrounds_averaged/#{semiuniprot}_#{species}.txt"
        next unless File.exist?(background_fn)

        output_folder = File.join("models/thresholds/#{model_type}/", "#{species}_background", uniprot)
        ext_glob = (model_type == 'mono') ? '*.pwm' : '*.dpwm'
        next  if Dir.exist?(output_folder) && FileList[File.join(pwm_folder, ext_glob)].pathmap('%n').sort == FileList[File.join(output_folder, '*.thr')].pathmap('%n').sort
        FileUtils.rm_rf output_folder
        FileUtils.mkdir_p output_folder
        Dir.glob("#{pwm_folder}/*").sort.each do |pwm_fn|
          output_fn = File.join(output_folder, File.basename(pwm_fn, File.extname(pwm_fn)) + '.thr')
          next  if File.exist?(output_fn)
          puts Ape.precalculate_thresholds_cmd pwm_fn,
                                               output_folder: output_folder,
                                               background: File.read(background_fn), # we always use dinucleotide background
                                               threshold_grid: ['1e-15', '1.0', '1.01', 'mul'],
                                               discretization: 1000,
                                               additional_options: additional_options,
                                               mode: :di
        end
      end
    end
  end
end
