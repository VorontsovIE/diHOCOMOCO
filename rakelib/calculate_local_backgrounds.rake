require 'fasta_sequence'
require 'sequence_dataset'
require 'frequencies'

desc 'Calculate local background for each dataset'
task :calculate_local_backgrounds do
  # mkdir_p 'control/local_backgrounds/mono'
  mkdir_p 'control/local_backgrounds/di'

  SequenceDataset.each_dataset {|dataset|
    # File.write(dataset.local_mono_background_path, dataset.local_mono_background)
    File.write(dataset.local_di_background_path, dataset.local_di_background)
  }
end

task :average_local_backgrounds do
  mkdir_p 'control/local_backgrounds_averaged/'
  Dir.glob('control/local_backgrounds/di/*.txt').group_by{|fn|
    File.basename(fn).split('.').first
  }.each{|uniprot, fns|
    bg = fns.map{|fn|
      File.read(fn).split(',').map(&:to_f)
    }.transpose.map{|dinuc_probs|
      dinuc_probs.inject(0.0, &:+) / dinuc_probs.size
    }
    File.write("control/local_backgrounds_averaged/#{uniprot}.txt", bg.join(','))
  }
end
