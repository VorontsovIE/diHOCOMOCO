require 'models'

module SequenceLogoGenerator
  def self.mono_cmd(pcm_file:, output_folder:, orientation: 'both', x_unit: 30, y_unit: 60, additional_options: [])
    opts = []
    opts += ['--logo-folder', output_folder]
    opts += ['--orientation', orientation]
    opts += ['--x-unit', x_unit.to_s]
    opts += ['--y-unit', y_unit.to_s]
    opts += additional_options.flatten.map(&:to_s)

    Shellwords.shelljoin(['sequence_logo', *opts, pcm_file])
  end

  def self.di_cmd(pcm_file:, output_folder:, orientation: 'both', x_unit: 30, y_unit: 60)
    ['direct', 'revcomp'].map do |orient|
      next  unless orientation == 'both' || orientation == orient
      basename = File.basename(pcm_file, File.extname(pcm_file))
      output_file = File.join(output_folder, "#{basename}_#{orient}.png")
      opts = [x_unit.to_s, y_unit.to_s]
      opts << '--revcomp'  if orient == 'revcomp'
      Shellwords.shelljoin(['ruby', './pmflogo/dpmflogo3.rb', pcm_file, output_file, *opts])
    end
  end

  def self.print_di_cmd_if_necessary(pcm_file:, output_folder:, orientation: 'both', x_unit: 30, y_unit: 60)
    motif = File.basename(pcm_file, File.extname(pcm_file))
    cmd = di_cmd( pcm_file: pcm_file, output_folder: output_folder,
                  orientation: orientation, x_unit: x_unit, y_unit: y_unit )
    if orientation == 'both'
      output_fn_1 = File.join(output_folder, "#{motif}_direct.png")
      output_fn_2 = File.join(output_folder, "#{motif}_revcomp.png")
      puts(cmd)  unless File.exist?(output_fn_1) && File.exist?(output_fn_2)
    else
      output_fn = File.join(output_folder, "#{motif}_#{orientation}.png")
      puts(cmd)  unless File.exist?(output_fn)
    end
  end

  def self.print_mono_cmd_if_necessary(pcm_file:, output_folder:, orientation: 'both', x_unit: 30, y_unit: 60, additional_options: [])
    motif = File.basename(pcm_file, File.extname(pcm_file))
    cmd = mono_cmd( pcm_file: pcm_file, output_folder: output_folder,
                    orientation: orientation, x_unit: x_unit, y_unit: y_unit,
                    additional_options: additional_options )
    if orientation == 'both'
      output_fn_1 = File.join(output_folder, "#{motif}_direct.png")
      output_fn_2 = File.join(output_folder, "#{motif}_revcomp.png")
      puts(cmd)  unless File.exist?(output_fn_1) && File.exist?(output_fn_2)
    else
      output_fn = File.join(output_folder, "#{motif}.png")
      puts(cmd)  unless File.exist?(output_fn)
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
  FileUtils.mkdir_p('final_collection/mono/logo/')
  FileUtils.mkdir_p('final_collection/mono/logo_large/')
  FileUtils.mkdir_p('final_collection/mono/logo_small/')
  Dir.glob('final_collection/mono/pcm/*.pcm').each do |fn|
    SequenceLogoGenerator.print_mono_cmd_if_necessary(
      pcm_file: fn,
      output_folder: 'final_collection/mono/logo/',
      **SequenceLogoGenerator::DefaultSizes[:big],
      additional_options: ['--no-threshold-lines']
    )

    SequenceLogoGenerator.print_mono_cmd_if_necessary(
      pcm_file: fn,
      output_folder: 'final_collection/mono/logo_large/',
      **SequenceLogoGenerator::DefaultSizes[:large],
      additional_options: ['--no-threshold-lines']
    )

    motif_length = ModelKind.get('mono').read_pcm(fn).length
    if motif_length > SequenceLogoGenerator::LongModelThreshold
      SequenceLogoGenerator.print_mono_cmd_if_necessary(
        pcm_file: fn,
        output_folder: 'final_collection/mono/logo_small/',
        **SequenceLogoGenerator::DefaultSizes[:small_for_long_models],
        additional_options: ['--no-threshold-lines']
      )
    else
      SequenceLogoGenerator.print_mono_cmd_if_necessary(
        pcm_file: fn,
        output_folder: 'final_collection/mono/logo_small/',
        **SequenceLogoGenerator::DefaultSizes[:small],
        additional_options: ['--no-threshold-lines']
      )
    end
  end
end


task "sequence_logos:di" do
  FileUtils.mkdir_p('final_collection/di/logo/')
  FileUtils.mkdir_p('final_collection/di/logo_large/')
  FileUtils.mkdir_p('final_collection/di/logo_small/')
  Dir.glob('final_collection/di/pcm/*.dpcm').each do |fn|
    SequenceLogoGenerator.print_di_cmd_if_necessary(
      pcm_file: fn,
      output_folder: 'final_collection/di/logo/',
      **SequenceLogoGenerator::DefaultSizes[:big]
    )

    SequenceLogoGenerator.print_di_cmd_if_necessary(
      pcm_file: fn,
      output_folder: 'final_collection/di/logo_large/',
      **SequenceLogoGenerator::DefaultSizes[:large]
    )


    motif_length = ModelKind.get('di').read_pcm(fn).length
    if motif_length > SequenceLogoGenerator::LongModelThreshold
      SequenceLogoGenerator.print_mono_cmd_if_necessary(
        pcm_file: fn,
        output_folder: 'final_collection/di/logo_small/',
        **SequenceLogoGenerator::DefaultSizes[:small_for_long_models],
        additional_options: ['--no-threshold-lines', '--dinucleotide']
      )
    else
      SequenceLogoGenerator.print_mono_cmd_if_necessary(
        pcm_file: fn,
        output_folder: 'final_collection/di/logo_small/',
        **SequenceLogoGenerator::DefaultSizes[:small],
        additional_options: ['--no-threshold-lines', '--dinucleotide']
      )
    end
  end
end
