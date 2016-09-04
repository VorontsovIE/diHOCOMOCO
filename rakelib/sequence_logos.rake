require 'models'

module SequenceLogoGenerator
  def self.run(pcm_files:, output_folder:, orientation: 'both', x_unit: 30, y_unit: 60, additional_options: [])
    return  if pcm_files.empty?

    FileUtils.mkdir_p(output_folder)  unless Dir.exist?(output_folder)

    opts = []
    opts += ['--logo-folder', output_folder]
    opts += ['--orientation', orientation]
    opts += ['--x-unit', x_unit.to_s]
    opts += ['--y-unit', y_unit.to_s]
    opts += additional_options.flatten.map(&:to_s)

    Open3.popen2('sequence_logo', *opts){|fread, fwrite|
      fread.puts Shellwords.join(pcm_files)
      fread.close
    }
  end

  def self.run_dinucleotide(pcm_files:, output_folder:, orientation: 'both', x_unit: 30, y_unit: 60)
    return  if pcm_files.empty?

    FileUtils.mkdir_p(output_folder)  unless Dir.exist?(output_folder)

    ['direct', 'revcomp'].each do |orient|
      next  unless orientation == 'both' || orientation == orient
      pcm_files.each do |pcm_file|
        basename = File.basename(pcm_file, File.extname(pcm_file))
        output_file = File.join(output_folder, "#{basename}_#{orient}.png")
        opts = [x_unit.to_s, y_unit.to_s]
        opts << '--revcomp'  if orient == 'revcomp'
        system 'ruby', './pmflogo/dpmflogo3.rb', pcm_file, output_file, *opts
      end
    end
  end


  DefaultSizes = {
    large: {x_unit: 100, y_unit: 200},
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
      output_folder: File.join('models/logo/', uniprot),
      **SequenceLogoGenerator::DefaultSizes[:big],
      additional_options: ['--no-threshold-lines']
    )

    SequenceLogoGenerator.run(
      pcm_files: models,
      output_folder: File.join('models/logo_large/', uniprot),
      **SequenceLogoGenerator::DefaultSizes[:large],
      additional_options: ['--no-threshold-lines']
    )

    SequenceLogoGenerator.run(
      pcm_files: normal_models,
      output_folder: File.join('models/logo_small/', uniprot),
      **SequenceLogoGenerator::DefaultSizes[:small],
      additional_options: ['--no-threshold-lines']
    )

    SequenceLogoGenerator.run(
      pcm_files: long_models,
      output_folder: File.join('models/logo_small/', uniprot),
      **SequenceLogoGenerator::DefaultSizes[:small_for_long_models],
      additional_options: ['--no-threshold-lines']
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

    SequenceLogoGenerator.run_dinucleotide(
      pcm_files: models,
      output_folder: File.join('models/logo/', uniprot),
      **SequenceLogoGenerator::DefaultSizes[:big]
    )

    SequenceLogoGenerator.run_dinucleotide(
      pcm_files: models,
      output_folder: File.join('models/logo_large/', uniprot),
      **SequenceLogoGenerator::DefaultSizes[:large]
    )

    SequenceLogoGenerator.run(
      pcm_files: normal_models,
      output_folder: File.join('models/logo_small/', uniprot),
      **SequenceLogoGenerator::DefaultSizes[:small],
      additional_options: ['--no-threshold-lines', '--dinucleotide']
    )

    SequenceLogoGenerator.run(
      pcm_files: long_models,
      output_folder: File.join('models/logo_small/', uniprot),
      **SequenceLogoGenerator::DefaultSizes[:small_for_long_models],
      additional_options: ['--no-threshold-lines', '--dinucleotide']
    )
  end
  task 'sequence_logos:di' => "sequence_logos:di:#{uniprot}"
end
