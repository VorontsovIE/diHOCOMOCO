require 'roc_curve'
require 'models'
require 'sequence_dataset'

desc 'Calculate model AUC on existing datasets'
task :calculate_auc => [:calculate_auc_mono, :calculate_auc_di]

desc 'Calculate model AUC on existing datasets for mononucleotide models'
task :calculate_auc_mono

desc 'Calculate model AUC on existing datasets for dinucleotide models'
task :calculate_auc_di


SequenceDataset.each_file_by_glob('control/control/*.mfa') do |control|
  task "calculate_auc_mono:#{control.name}" do
    
    output_dir = File.join('occurences/auc/mono/', control.uniprot) 
    mkdir_p(output_dir)  unless Dir.exist?(output_dir)
    output_fn = File.join(output_dir, "#{control.name}.txt")
    
    File.open(output_fn, 'w') do |fw|
      Models.mono_models_by_uniprot(control.uniprot).each do |model|
        corrected_pvalues_fn = File.join('occurences/corrected_pvalues/mono/', control.uniprot, model.full_name, "#{control.name}.txt")
        corrected_pvalues = File.readlines(corrected_pvalues_fn).map(&:to_f)
        roc = roc_curve(corrected_pvalues)
        auc = calculate_auc(roc)
        fw.puts [model.full_name, auc].join("\t")
      end
    end
  end
  task :calculate_auc_mono => "calculate_auc_mono:#{control.name}"


  task "calculate_auc_di:#{control.name}" do
    
    output_dir = File.join('occurences/auc/di/', control.uniprot) 
    mkdir_p(output_dir)  unless Dir.exist?(output_dir)
    output_fn = File.join(output_dir, "#{control.name}.txt")
    
    File.open(output_fn, 'w') do |fw|
      Models.di_models_by_uniprot(control.uniprot).each do |model|
        corrected_pvalues_fn = File.join('occurences/corrected_pvalues/di/', control.uniprot, model.full_name, "#{control.name}.txt")
        corrected_pvalues = File.readlines(corrected_pvalues_fn).map(&:to_f)
        roc = roc_curve(corrected_pvalues)
        auc = calculate_auc(roc)
        fw.puts [model.full_name, auc].join("\t")
      end
    end
  end
  task :calculate_auc_di => "calculate_auc_di:#{control.name}"
end
