require 'rake/file_utils_ext'

module Ape
  def self.precalculate_thresholds_cmd(models_folder,
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
      
    FileUtils.mkdir_p output_folder  unless Dir.exist?(output_folder)

    cmd = ['java', '-Xmx1G', '-cp', 'ape.jar', package]

    opts = []
    opts += ['--silent']
    opts += ['--pvalues', threshold_grid.join(',')]  if threshold_grid
    opts += ['--discretization', discretization.to_s]  if discretization
    opts += ['--background', background.to_s]  if background
    opts += additional_options
    
    [*cmd, models_folder, output_folder, *opts].shelljoin
  end
end
