require 'models'

desc 'Draw sequence logos for all motifs'
task  :sequence_logos do
  Models.mono_uniprots.each do |uniprot|
    pcm_files = Dir.glob(File.join('models/pcm/mono/all/', uniprot, '*.pcm'))

    next  if pcm_files.empty?
    $stderr.puts pcm_files

    logo_folder = File.join('models/logo/', uniprot)
    mkdir_p logo_folder  unless Dir.exist?(logo_folder)

    opts = []
    opts += ['--logo-folder', logo_folder]
    Open3.popen2('sequence_logo', *opts){|fread, fwrite|
      fread.puts pcm_files
      fread.close
    }
  end


  Models.di_uniprots.each do |uniprot|
    pcm_files = Dir.glob(File.join('models/pcm/di/all/', uniprot, '*.dpcm'))

    next  if pcm_files.empty?
    $stderr.puts pcm_files

    logo_folder = File.join('models/logo/', uniprot)
    mkdir_p logo_folder  unless Dir.exist?(logo_folder)

    opts = []
    opts += ['--logo-folder', logo_folder]
    opts += ['--dinucleotide']
    Open3.popen2('sequence_logo', *opts){|fread, fwrite|
      fread.puts pcm_files
      fread.close
    }
  end
end
