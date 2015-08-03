module Ape
  def self.run_precalculate_thresholds(models_folder, 
                                      output_folder:,
                                      background: nil,
                                      threshold_grid: nil,
                                      discretization: nil,
                                      additional_options: [],
                                      mode:)
    case mode
    when /^mono$/
      package = 'ru.autosome.ape.PrecalculateThresholds'
    when /^di$/
      package = 'ru.autosome.ape.di.PrecalculateThresholds'
    else
      raise "Unknown mode `#{mode}`"
    end
      
    mkdir_p output_folder  unless Dir.exist?(output_folder)

    cmd = ['java', '-Xmx1G', '-cp', 'ape-2.0.1.jar', package]
    opts = []
    opts += ['--silent']
    opts += ['--pvalues', threshold_grid.join(',')]  if threshold_grid
    opts += ['--discretization', discretization.to_s]  if discretization
    opts += ['--background', background.to_s]  if background
    opts += additional_options
    
    sh *cmd, control_filename, model_filename, *opts, {out: output_file}, {}
    sh *script_cmd, motifs_folder, output_folder, *opts
  end
end
