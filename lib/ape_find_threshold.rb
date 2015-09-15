require 'open3'

module Ape
  def self.run_find_threshold(model_filename, pvalues,
                          background: nil,
                          discretization: nil,
                          additional_options: [],
                          mode:)
    case mode
    when /^mono$/
      package = 'ru.autosome.ape.FindThreshold'
    when /^di$/
      package = 'ru.autosome.ape.di.FindThreshold'
    else
      raise "Unknown mode `#{mode}`"
    end

    cmd = ['java', '-Xmx1G', '-cp', 'ape.jar', package]
    opts = []
    opts += ['--discretization', discretization.to_s]  if discretization
    opts += ['--background', background.to_s]  if background
    opts += additional_options

    result = nil
    $stderr.puts [*cmd, model_filename, *opts].join(' ')
    Open3.popen2(*cmd, model_filename, *opts){|fread, fwrite|
      fread.puts pvalues
      fread.close
      result = fwrite.read
    }

    result.lines.reject{|line|
      line.start_with?('#')
    }.map{|line|
      line.chomp.split("\t")
    }.map{|row|
      [row[0].to_f, row[3].to_f]
    }.to_h
  end
end
