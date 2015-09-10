require 'models'

desc 'Draw sequence logos for all motifs'
task :sequence_logos => ['sequence_logos:mono', 'sequence_logos:di']

Models.mono_uniprots.each do |uniprot|
  task "sequence_logos:mono:#{uniprot}" do
    pcm_files = FileList[File.join('models/pcm/mono/all/', uniprot, '*.pcm')]
    next  if pcm_files.empty?

    logo_folder = File.join((ENV['logo_folder'] || 'models/logo/'), uniprot)
    mkdir_p logo_folder  unless Dir.exist?(logo_folder)

    opts = []
    opts += ['--logo-folder', logo_folder]
    opts += ['--orientation', 'both']

    x_unit = (ENV['x_unit'] || 30).to_i
    y_unit = (ENV['y_unit'] || 60).to_i
    opts += ['--x-unit', x_unit.to_s]
    opts += ['--y-unit', y_unit.to_s]
    Open3.popen2('sequence_logo', *opts){|fread, fwrite|
      fread.puts Shellwords.join(pcm_files)
      fread.close
    }
  end
  task 'sequence_logos:mono' => "sequence_logos:mono:#{uniprot}"
end


Models.di_uniprots.each do |uniprot|
  task "sequence_logos:di:#{uniprot}" do
    pcm_files = FileList[File.join('models/pcm/di/all/', uniprot, '*.dpcm')]
    next  if pcm_files.empty?

    logo_folder = File.join((ENV['logo_folder'] || 'models/logo/'), uniprot)
    mkdir_p logo_folder  unless Dir.exist?(logo_folder)

    opts = []
    opts += ['--logo-folder', logo_folder]
    opts += ['--dinucleotide']
    opts += ['--orientation', 'both']

    x_unit = (ENV['x_unit'] || 30).to_i
    y_unit = (ENV['y_unit'] || 60).to_i
    opts += ['--x-unit', x_unit.to_s]
    opts += ['--y-unit', y_unit.to_s]
    Open3.popen2('sequence_logo', *opts){|fread, fwrite|
      fread.puts Shellwords.join(pcm_files)
      fread.close
    }
  end
  task 'sequence_logos:di' => "sequence_logos:di:#{uniprot}"
end
