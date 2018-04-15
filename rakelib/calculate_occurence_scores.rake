require 'models'
require 'shellwords'

# Takes sequences from fasta files and rename them so that each sequence is named in the form {DatasetName}:{SequenceLength}
def join_fasta_renamed(mfa_files, output_stream:)
  mfa_files.each{|fn|
    control_name = File.basename(fn, '.mfa')
    File.readlines(fn).each_slice(2){|_hdr, seq|
      output_stream.puts ">#{control_name}:#{seq.strip.length}\n#{seq}"
    }
  }
end

desc 'Calculate scores for each model on each control'
task :calculate_occurence_scores => [:calculate_occurence_scores_mono, :calculate_occurence_scores_di]

desc 'Calculate scores for each mononucleotide model on each control'
task :calculate_occurence_scores_mono

desc 'Calculate scores for each dinucleotide model on each control'
task :calculate_occurence_scores_di

desc 'group controls for a protein into single file'
task "group_controls" do
  FileUtils.mkdir_p "control/grouped"
  factors = Models.all_models.map(&:semiuniprot).uniq.sort
  ['HUMAN', 'MOUSE'].each do |species| # dataset characteristic
    factors.each do |factor|
      result_fn = "control/grouped/#{factor}_#{species}.mfa"
      tf_controls = Dir.glob("control/control/#{factor}_#{species}.*.mfa")
      next  if File.exist?(result_fn)
      next  if tf_controls.empty?
      File.open(result_fn, 'w'){|fw|  join_fasta_renamed(tf_controls, output_stream: fw)  }
    end
  end
end

[ModelKind::Mono.new, ModelKind::Di.new].each do |model_kind|
  task "calculate_occurence_scores_#{model_kind}" do
    models_by_factor = Models.models_by_type(model_kind).group_by(&:semiuniprot)
    ['HUMAN', 'MOUSE'].each do |species| # dataset characteristic; models can be of any species
      FileUtils.mkdir_p "auc/#{model_kind}/#{species}_datasets"

      models_by_factor.each{|factor, factor_models|
        factor_models.each{|model|
          control_fn = "control/grouped/#{factor}_#{species}.mfa"
          thresholds_fn = "models/thresholds/#{model_kind}/#{species}_background/#{model.uniprot}/#{model.full_name}.thr"
          output_fn = "auc/#{model_kind}/#{species}_datasets/#{model.full_name}.txt"
          next  unless File.exist?(thresholds_fn)
          next  if File.exist?(output_fn)
          puts "java -cp sarus.jar #{model_kind.sarus_class} #{control_fn} #{model.path_to_pwm} besthit  --output-scoring-mode pvalue  --pvalues-file #{thresholds_fn}" + \
               " | ruby calculate_auc.rb #{model.length}" + \
               " > #{output_fn}"
        }
      }
    end
  end
end
