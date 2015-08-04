require 'stringio'
require 'rake/file_utils_ext'

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
    
    scores_io = StringIO.new(scores.map(&:to_s).join("\n"))
    Rake::FileUtilsExt.sh *cmd, model_filename, *opts, {out: output_file, in: scores_io}, {}
  end
end
