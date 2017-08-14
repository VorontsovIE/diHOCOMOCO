require 'models'

module SequenceLogoGenerator
  def self.mono_cmd(pcm_file:, output_folder:, orientation: 'both', x_unit: 30, y_unit: 60, additional_options: [])
    FileUtils.mkdir_p(output_folder)  unless Dir.exist?(output_folder)

    opts = []
    opts += ['--logo-folder', output_folder]
    opts += ['--orientation', orientation]
    opts += ['--x-unit', x_unit.to_s]
    opts += ['--y-unit', y_unit.to_s]
    opts += additional_options.flatten.map(&:to_s)

    Shellwords.shelljoin(['sequence_logo', *opts, pcm_file])
  end

  def self.di_cmds(pcm_file:, output_folder:, orientation: 'both', x_unit: 30, y_unit: 60)
    FileUtils.mkdir_p(output_folder)  unless Dir.exist?(output_folder)

    ['direct', 'revcomp'].map do |orient|
      next  unless orientation == 'both' || orientation == orient
      basename = File.basename(pcm_file, File.extname(pcm_file))
      output_file = File.join(output_folder, "#{basename}_#{orient}.png")
      opts = [x_unit.to_s, y_unit.to_s]
      opts << '--revcomp'  if orient == 'revcomp'
      Shellwords.shelljoin(['ruby', './pmflogo/dpmflogo3.rb', pcm_file, output_file, *opts])
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

task "sequence_logos:mono" do
  Dir.glob('final_collection/mono/pcm/*.pcm').each do |fn|
    puts SequenceLogoGenerator.mono_cmd(
      pcm_file: fn,
      output_folder: 'final_collection/mono/logo/',
      **SequenceLogoGenerator::DefaultSizes[:big],
      additional_options: ['--no-threshold-lines']
    )

    puts SequenceLogoGenerator.mono_cmd(
      pcm_file: fn,
      output_folder: 'final_collection/mono/logo_large/',
      **SequenceLogoGenerator::DefaultSizes[:large],
      additional_options: ['--no-threshold-lines']
    )

    motif_length = ModelKind.get('mono').read_pcm(fn).length
    if motif_length > SequenceLogoGenerator::LongModelThreshold
      puts SequenceLogoGenerator.mono_cmd(
        pcm_file: fn,
        output_folder: 'final_collection/mono/logo_small/',
        **SequenceLogoGenerator::DefaultSizes[:small_for_long_models],
        additional_options: ['--no-threshold-lines']
      )
    else
      puts SequenceLogoGenerator.mono_cmd(
        pcm_file: fn,
        output_folder: 'final_collection/mono/logo_small/',
        **SequenceLogoGenerator::DefaultSizes[:small],
        additional_options: ['--no-threshold-lines']
      )
    end
  end
end


task "sequence_logos:di" do
  Dir.glob('final_collection/di/pcm/*.dpcm').each do |fn|
    puts SequenceLogoGenerator.di_cmds(
      pcm_file: fn,
      output_folder: 'final_collection/di/logo/',
      **SequenceLogoGenerator::DefaultSizes[:big]
    )

    puts SequenceLogoGenerator.di_cmds(
      pcm_file: fn,
      output_folder: 'final_collection/di/logo_large/',
      **SequenceLogoGenerator::DefaultSizes[:large]
    )

    motif_length = ModelKind.get('di').read_pcm(fn).length
    if motif_length > SequenceLogoGenerator::LongModelThreshold
      puts SequenceLogoGenerator.mono_cmd(
        pcm_file: fn,
        output_folder: 'final_collection/di/logo_small/',
        **SequenceLogoGenerator::DefaultSizes[:small_for_long_models],
        additional_options: ['--no-threshold-lines', '--dinucleotide']
      )
    else
      puts SequenceLogoGenerator.mono_cmd(
        pcm_file: fn,
        output_folder: 'final_collection/di/logo_small/',
        **SequenceLogoGenerator::DefaultSizes[:small],
        additional_options: ['--no-threshold-lines', '--dinucleotide']
      )
    end
  end
end
