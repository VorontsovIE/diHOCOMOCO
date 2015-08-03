# for `besthit suppress` mode
SarusResults = Struct.new(:score, :position, :orientation) do
  def self.from_string(str)
    score, pos, orientation = str.chomp.split("\t")
    self.new(score.to_f, pos.to_i, orientation.to_sym)
  end

  def self.each_in_file(filename)
    return enum_for(:each_in_file, filename)  unless block_given?
    File.open(filename) do |f|
      f.each_line{|line| yield self.from_string(line) }
    end
  end
end

module Sarus
  def self.run_besthits(control_filename, model_filename, output_file:, mode:)
    case mode
    when /^mono$/
      package = 'ru.autosome.SARUS'
    when /^di$/
      package = 'ru.autosome.di.SARUS'
    else
      raise "Unknown mode `#{mode}`"
    end
      
    output_folder = File.dirname(output_file)
    mkdir_p output_folder  unless Dir.exist?(output_folder)

    cmd = ['java', '-Xmx1G', '-cp', 'sarus.jar', package]
    opts = ['besthit', 'suppress']
    sh *cmd, control_filename, model_filename, *opts, {out: output_file}, {}
  end
end
