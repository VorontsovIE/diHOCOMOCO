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

  DefaultSizes = {
    big: {x_unit: 30, y_unit: 60},
    small: {x_unit: 15, y_unit: 30},
    small_for_long_models: {x_unit: 10, y_unit: 20},
  }
  LongModelThreshold = 20
end

desc 'Draw sequence logos for all motifs'
task :sequence_logos => ['sequence_logos:mono', 'sequence_logos:di']

Models.mono_uniprots.each do |uniprot|
  task "sequence_logos:mono:#{uniprot}" do
    models = FileList[File.join('models/pcm/mono/all/', uniprot, '*.pcm')]
    normal_models = models.select{|fn|
      ModelKind.get('mono').read_pcm(fn).length <= SequenceLogoGenerator::LongModelThreshold
    }
    long_models = models.select{|fn|
      ModelKind.get('mono').read_pcm(fn).length > SequenceLogoGenerator::LongModelThreshold
    }

    SequenceLogoGenerator.run(
      pcm_files: models,
      output_folder: 'models/logo/',
      **SequenceLogoGenerator::DefaultSizes[:big]
    )

    SequenceLogoGenerator.run(
      pcm_files: normal_models,
      output_folder: 'models/logo_small/',
      **SequenceLogoGenerator::DefaultSizes[:small]
    )

    SequenceLogoGenerator.run(
      pcm_files: long_models,
      output_folder: 'models/logo_small/',
      **SequenceLogoGenerator::DefaultSizes[:small_for_long_models]
    )
  end
  task 'sequence_logos:mono' => "sequence_logos:mono:#{uniprot}"
end


Models.di_uniprots.each do |uniprot|
  task "sequence_logos:di:#{uniprot}" do
    models = FileList[File.join('models/pcm/di/all/', uniprot, '*.dpcm')]
    normal_models = models.select{|fn|
      ModelKind.get('di').read_pcm(fn).length <= SequenceLogoGenerator::LongModelThreshold
    }
    long_models = models.select{|fn|
      ModelKind.get('di').read_pcm(fn).length > SequenceLogoGenerator::LongModelThreshold
    }

    SequenceLogoGenerator.run(
      pcm_files: models,
      output_folder: 'models/logo/',
      **SequenceLogoGenerator::DefaultSizes[:big]
      additional_options: ['--dinucleotide']
    )

    SequenceLogoGenerator.run(
      pcm_files: normal_models,
      output_folder: 'models/logo_small/',
      **SequenceLogoGenerator::DefaultSizes[:small]
      additional_options: ['--dinucleotide']
    )

    SequenceLogoGenerator.run(
      pcm_files: long_models,
      output_folder: 'models/logo_small/',
      **SequenceLogoGenerator::DefaultSizes[:small_for_long_models]
      additional_options: ['--dinucleotide']
    )
  end
  task 'sequence_logos:di' => "sequence_logos:di:#{uniprot}"
end
