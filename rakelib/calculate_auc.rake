require 'roc_curve'
require 'models'
require 'sequence_dataset'

desc 'Calculate model AUC on existing datasets'
task :calculate_auc => [:calculate_auc_mono, :calculate_auc_di]

task :calculate_auc_mono do
  uniprots = SequenceDataset.each_file_by_glob('control/control/*.mfa').map(&:uniprot).uniq.sort

  uniprots.each do |uniprot|
    output_dir = File.join('occurences/auc/mono/', uniprot)
    mkdir_p(output_dir)  unless Dir.exist?(output_dir)

    Models.mono_models_by_uniprot(uniprot).each do |model|
      output_fn = File.join(output_dir, "#{model.full_name}.txt")

      File.open(output_fn, 'w') do |fw|
        SequenceDataset.each_file_by_glob("control/control/#{uniprot}^*.mfa") do |control|
          corrected_pvalues_fn = File.join('occurences/corrected_pvalues/mono/', control.uniprot, model.full_name, "#{control.name}.txt")
          corrected_pvalues = File.readlines(corrected_pvalues_fn).map(&:to_f)
          roc = roc_curve(corrected_pvalues)
          auc = calculate_auc(roc)
          fw.puts [control.name, auc].join("\t")
        end
      end

    end
  end
end

task :calculate_auc_di do
  uniprots = SequenceDataset.each_file_by_glob('control/control/*.mfa').map(&:uniprot).uniq.sort

  uniprots.each do |uniprot|
    output_dir = File.join('occurences/auc/di/', uniprot)
    mkdir_p(output_dir)  unless Dir.exist?(output_dir)

    Models.di_models_by_uniprot(uniprot).each do |model|
      output_fn = File.join(output_dir, "#{model.full_name}.txt")

      File.open(output_fn, 'w') do |fw|
        SequenceDataset.each_file_by_glob("control/control/#{uniprot}^*.mfa") do |control|
          corrected_pvalues_fn = File.join('occurences/corrected_pvalues/di/', control.uniprot, model.full_name, "#{control.name}.txt")
          corrected_pvalues = File.readlines(corrected_pvalues_fn).map(&:to_f)
          roc = roc_curve(corrected_pvalues)
          auc = calculate_auc(roc)
          fw.puts [control.name, auc].join("\t")
        end
      end

    end
  end
end
