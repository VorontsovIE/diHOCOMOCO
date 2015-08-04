require 'rake/file_utils_ext'
require 'open3'

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

    output_folder = File.dirname(output_file)
    FileUtils.mkdir_p output_folder  unless Dir.exist?(output_folder)

    cmd = ['java', '-Xmx1G', '-cp', 'ape.jar', package]
    opts = []
    opts += ['--discretization', discretization.to_s]  if discretization
    opts += ['--background', background.to_s]  if background
    opts += additional_options
    
    $stderr.puts [*cmd, model_filename, *opts].join(' ')
    Open3.popen2(*cmd, model_filename, *opts){|fread, fwrite|
      fread.puts scores
      fread.close
      File.write(output_file, fwrite.read)
    }
  end
end
