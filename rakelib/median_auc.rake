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

  output_dir = File.join('occurences/median_auc/')
  mkdir_p(output_dir)  unless Dir.exist?(output_dir)

  SequenceDataset.each_uniprot do |uniprot|
    median_auc_by_model = Models.all_models_by_uniprot(uniprot).map{|model|
      model_auc_fn = File.join('occurences/auc/', uniprot, "#{model.full_name}.txt")
      [model, median_auc_in_file(model_auc_fn)]
    }.to_h
    print_sorted_hash_to_file(File.join(output_dir, "#{uniprot}.txt"), median_auc_by_model)
  end
end
