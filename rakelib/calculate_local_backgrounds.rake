require 'weighted_sequence'
require 'sequence_dataset'
require 'background'

desc 'Calculate local background for each dataset'
task :calculate_local_backgrounds do
  mkdir_p 'control/local_backgrounds/mono'
  mkdir_p 'control/local_backgrounds/di'

  SequenceDataset.each_by_glob('control/control/*.mfa') {|dataset|
    File.write("control/local_backgrounds/mono/#{dataset.name}.txt", dataset.local_mono_background)
    File.write("control/local_backgrounds/di/#{dataset.name}.txt", dataset.local_di_background)
  }
end
