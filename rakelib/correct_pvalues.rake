require 'models'
require 'sarus_results'
require 'ape_find_pvalue'
require 'find_pvalue_results'

# from_file - is a threshold-pvalue list
# to_file - is a list of corrected pvalues
def correct_pvalues(from_file, to_file, median_length:, model_length:)
  pvalues = FindPvalueResults.each_in_file(from_file).map(&:pvalue)
  corrected_pvalues = pvalues.map{|pvalue|
    1.0 - (1.0 - pvalue) ** (2 * (median_length - model_length + 1))
  }
  output_folder = File.dirname(to_file)
  mkdir_p(output_folder)  unless Dir.exist?(output_folder)
  File.write(to_file, corrected_pvalues.join("\n"))
end

desc 'Correct pvalues for sequence length'
task :correct_pvalues => [:corrected_pvalues_mono, :corrected_pvalues_di]

desc 'Correct pvalues for sequence length for mononucleotide models'
task :corrected_pvalues_mono

desc 'Correct pvalues for sequence length for dinucleotide models'
task :corrected_pvalues_di


SequenceDataset.each_dataset do |control|
  task "corrected_pvalues_mono:#{control.name}" do
    lengths = control.each_sequence.map(&:length)
    median_length = lengths.sort[lengths.size / 2]
    Models.mono_models_by_uniprot(control.uniprot).each do |model|
      pvalues_fn = File.join('occurences/pvalues/mono/', control.uniprot, model.full_name, "#{control.name}.txt")
      corrected_pvalues_fn = File.join('occurences/corrected_pvalues/mono/', control.uniprot, model.full_name, "#{control.name}.txt")
      correct_pvalues(pvalues_fn, corrected_pvalues_fn, median_length: median_length, model_length: model.length)
    end
  end
  task :corrected_pvalues_mono => "corrected_pvalues_mono:#{control.name}"


  task "corrected_pvalues_di:#{control.name}" do
    lengths = control.each_sequence.map(&:length)
    median_length = lengths.sort[lengths.size / 2]
    Models.di_models_by_uniprot(control.uniprot).each do |model|
      pvalues_fn = File.join('occurences/pvalues/di/', control.uniprot, model.full_name, "#{control.name}.txt")
      corrected_pvalues_fn = File.join('occurences/corrected_pvalues/di/', control.uniprot, model.full_name, "#{control.name}.txt")
      correct_pvalues(pvalues_fn, corrected_pvalues_fn, median_length: median_length, model_length: model.length)
    end
  end
  task :corrected_pvalues_di => "corrected_pvalues_di:#{control.name}"
end
