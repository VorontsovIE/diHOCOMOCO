require 'median'
require 'models'
require 'sequence_dataset'

def median_auc_in_file(filename)
  model_aucs = File.readlines(filename).map{|line|
    control_name, auc = line.chomp.split("\t")
    auc.to_f
  }
  median(model_aucs)
end

def print_sorted_hash_to_file(filename, hsh)
  File.open(filename, 'w') do |fw|
    hsh.sort_by{|model, median_auc|
      median_auc
    }.reverse_each{|model, median_auc|
      fw.puts [model.full_name, median_auc].join("\t")
    }
  end
end

desc 'Calculate model AUC on existing datasets'
task :median_auc do

  output_dir_mono = File.join('occurences/median_auc/mono/')
  mkdir_p(output_dir_mono)  unless Dir.exist?(output_dir_mono)

  output_dir_di = File.join('occurences/median_auc/di/')
  mkdir_p(output_dir_di)  unless Dir.exist?(output_dir_di)

  SequenceDataset.each_uniprot do |uniprot|
    median_auc_by_model_mono = Models.mono_models_by_uniprot(uniprot).map{|model|
      model_auc_fn = File.join('occurences/auc/mono/', uniprot, "#{model.full_name}.txt")
      [model, median_auc_in_file(model_auc_fn)]
    }.to_h

    print_sorted_hash_to_file(File.join(output_dir_mono, "#{uniprot}.txt"), median_auc_by_model_mono)

    median_auc_by_model_di = Models.di_models_by_uniprot(uniprot).map{|model|
      model_auc_fn = File.join('occurences/auc/di/', uniprot, "#{model.full_name}.txt")
      [model, median_auc_in_file(model_auc_fn)]
    }.to_h
    print_sorted_hash_to_file(File.join(output_dir_di, "#{uniprot}.txt"), median_auc_by_model_di)
  end
end
