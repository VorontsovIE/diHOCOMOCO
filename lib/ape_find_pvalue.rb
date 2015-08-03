module Ape
  def self.run_find_pvalue(model_filename, scores,
                          output_file:,
                          background: nil,
                          discretization: nil,
                          additional_options: [],
                          mode:)
    case mode
    when /^mono$/
      package = 'ru.autosome.ape.FindPvalue'
    when /^di$/
      package = 'ru.autosome.ape.di.FindPvalue'
    else
      raise "Unknown mode `#{mode}`"
    end
      
    mkdir_p output_folder  unless Dir.exist?(output_folder)

    cmd = ['java', '-Xmx1G', '-cp', 'ape-2.0.1.jar', package]
    opts = []
    opts += ['--discretization', discretization.to_s]  if discretization
    opts += ['--background', background.to_s]  if background
    opts += additional_options
    
    sh *cmd, model_filename, *scores.map(&:to_s), *opts, {out: output_file}, {}
  end
end
