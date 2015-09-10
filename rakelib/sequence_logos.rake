require 'models'

module SequenceLogoGenerator
  def self.run(pcm_files:, output_folder:, orientation: 'both', x_unit: 30, y_unit: 60, additional_options: [])
    return  if pcm_files.empty?

    output_subfolder = File.join(output_folder, uniprot)
    FileUtils.mkdir_p(output_subfolder)  unless Dir.exist?(output_subfolder)

    opts = []
    opts += ['--logo-folder', output_subfolder]
    opts += ['--orientation', orientation]
    opts += ['--x-unit', x_unit.to_s]
    opts += ['--y-unit', y_unit.to_s]
    opts += additional_options.flatten.map(&:to_s)

    Open3.popen2('sequence_logo', *opts){|fread, fwrite|
      fread.puts Shellwords.join(pcm_files)
      fread.close
    }
  end

end

desc 'Draw sequence logos for all motifs'
task :sequence_logos => ['sequence_logos:mono', 'sequence_logos:di']

Models.mono_uniprots.each do |uniprot|
  task "sequence_logos:mono:#{uniprot}" do
    SequenceLogoGenerator.run(
      pcm_files: FileList[File.join('models/pcm/mono/all/', uniprot, '*.pcm')],
      output_folder: 'models/logo/',
    )

    SequenceLogoGenerator.run(
      pcm_files: FileList[File.join('models/pcm/mono/all/', uniprot, '*.pcm')],
      output_folder: 'models/logo_small/',
      x_unit: 15, y_unit: 30,
    )
  end
  task 'sequence_logos:mono' => "sequence_logos:mono:#{uniprot}"
end


Models.di_uniprots.each do |uniprot|
  task "sequence_logos:di:#{uniprot}" do
    SequenceLogoGenerator.run(
      pcm_files: FileList[File.join('models/pcm/di/all/', uniprot, '*.dpcm')],
      output_folder: 'models/logo/',
      x_unit: 30, y_unit: 60,
      additional_options: ['--dinucleotide']
    )
    SequenceLogoGenerator.run(
      pcm_files: FileList[File.join('models/pcm/di/all/', uniprot, '*.dpcm')],
      output_folder: 'models/logo_small/',
      x_unit: 15, y_unit: 30,
      additional_options: ['--dinucleotide']
    )
  end
  task 'sequence_logos:di' => "sequence_logos:di:#{uniprot}"
end
