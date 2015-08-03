require 'weighted_sequence'
require 'sequence_dataset'
require 'frequencies'

desc 'Calculate local background for each dataset'
task :calculate_local_backgrounds do
  # mkdir_p 'control/local_backgrounds/mono'
  mkdir_p 'control/local_backgrounds/di'

  SequenceDataset.each_file_by_glob('control/control/*.mfa') {|dataset|
    # File.write(dataset.local_mono_background_path, dataset.local_mono_background)
    File.write(dataset.local_di_background_path, dataset.local_di_background)
  }
end
