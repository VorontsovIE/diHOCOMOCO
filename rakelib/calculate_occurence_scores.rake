require 'models'

desc 'Calculate scores for each model on each control'
task :calculate_occurence_scores => [:calculate_occurence_scores_mono, :calculate_occurence_scores_di]

desc 'Calculate scores for each mononucleotide model on each control'
task :calculate_occurence_scores_mono

desc 'Calculate scores for each dinucleotide model on each control'
task :calculate_occurence_scores_di

['mono', 'di'].each do |model_type|
  task "calculate_occurence_scores_#{model_type}" do
    sarus_class = (model_type == 'mono') ? 'ru.autosome.SARUS' : 'ru.autosome.di.SARUS'
    ['HUMAN', 'MOUSE'].each do |species| # dataset characteristic; models can be of any species
      FileUtils.mkdir_p "auc/#{model_type}/#{species}_datasets"
      Dir.glob("models/pwm/#{model_type}/all/*/*.{pwm,dpwm}").map{|fn|
        File.basename(fn, File.extname(fn))
      }.uniq.sort.map{|model_name|
        Model.new(model_name, model_type.to_sym)
      }.each do |model|
        thresholds_fn = "models/thresholds/#{model_type}/#{species}_background/#{model.uniprot}/#{model.full_name}.thr"
        output_fn = "auc/#{model_type}/#{species}_datasets/#{model.full_name}.txt"
        factor = model.uniprot.split('_').first # species-independent
        next  unless File.exist?(thresholds_fn)
        next  if File.exist?(output_fn)
        puts "find control/control/ -name '#{factor}_#{species}.*.mfa' " + \
             " | xargs -r -n1 -I{} basename -s .mfa '{}' " + \
             " | xargs -r -n1 -I{} echo \"cat control/control/{}.mfa | ruby renameMultifastaSequences.rb '{}'\" " + \
             " | bash " + \
             " | java -cp sarus.jar #{sarus_class} - #{model.path_to_pwm} besthit " + \
             " | ruby calculate_auc.rb #{model.length} #{thresholds_fn} #{model_type} " + \
             " > #{output_fn}"
      end
    end
  end
end
