require 'median'
require 'models'
require 'sarus_results'
require 'ape_find_pvalue'
require 'find_pvalue_results'

# from_file - is a threshold-pvalue list
# to_file - is a list of corrected pvalues
def correct_pvalues(from_file, to_file, median_length:, model_length:)
  return  if File.exist?(to_file)
  $stderr.puts "Correct pvalues: #{from_file} --> #{to_file}"
  pvalues = FindPvalueResults.each_in_file(from_file).map(&:pvalue)
  corrected_pvalues = pvalues.map{|pvalue|
    1.0 - (1.0 - pvalue) ** (2 * (median_length - model_length + 1))
  }
  output_folder = File.dirname(to_file)
  mkdir_p(output_folder)  unless Dir.exist?(output_folder)
  File.write(to_file, corrected_pvalues.join("\n"))
end

desc 'Correct pvalues for sequence length'
task :correct_pvalues do
  SequenceDataset.each_dataset do |control|
    lengths = control.each_sequence.map(&:length)
    Models.all_models_by_uniprot(control.uniprot).each do |model|
      correct_pvalues(File.join('occurences/pvalues/', control.uniprot, model.full_name, "#{control.name}.txt"),
                      File.join('occurences/corrected_pvalues/', control.uniprot, model.full_name, "#{control.name}.txt"),
                      median_length: median(lengths),
                      model_length: model.length)
    end
  end
end
