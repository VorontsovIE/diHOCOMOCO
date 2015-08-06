require 'roc_curve'
require 'models'
require 'sequence_dataset'

desc 'Calculate model AUC on existing datasets'
task :calculate_auc do
  SequenceDataset.each_uniprot do |uniprot|
    output_dir = File.join('occurences/auc/', uniprot)
    mkdir_p(output_dir)  unless Dir.exist?(output_dir)

    Models.all_models_by_uniprot(uniprot).each do |model|
      output_fn = File.join(output_dir, "#{model.full_name}.txt")

      File.open(output_fn, 'w') do |fw|
        SequenceDataset.each_for_uniprot(uniprot) do |control|
          corrected_pvalues_fn = File.join('occurences/corrected_pvalues/', control.uniprot, model.full_name, "#{control.name}.txt")
          corrected_pvalues = File.readlines(corrected_pvalues_fn).map(&:to_f)
          roc = roc_curve(corrected_pvalues)
          auc = calculate_auc(roc)
          fw.puts [control.name, auc].join("\t")
        end
      end

    end
  end
end
