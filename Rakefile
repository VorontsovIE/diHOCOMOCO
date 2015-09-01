$:.unshift File.absolute_path('lib')
require 'rake/clean'

desc 'Move all motifs/words into their folders, normalize their names, convert to PWMs and so on'
task 'collect_and_normalize_data' do
  Rake::Task['collect_and_normalize_data:collect_pcm'].invoke
  Rake::Task['collect_and_normalize_data:rename_motifs'].invoke
  Rake::Task['collect_and_normalize_data:convert_pcm_to_pwm'].invoke
end

task :default do
  Rake::Task['unpack'].invoke
  Rake::Task['collect_and_normalize_data'].invoke
  Rake::Task['remove_non_currated_models'].invoke

  Rake::Task['calculate_local_backgrounds'].invoke
  Rake::Task['precalculate_thresholds'].invoke

  Rake::Task['calculate_occurence_scores'].invoke
  Rake::Task['scores_to_pvalues'].invoke
  Rake::Task['correct_pvalues'].invoke
  Rake::Task['calculate_auc'].invoke

  Rake::Task['filter_motifs_and_datasets'].invoke
  Rake::Task['sequence_logos'].invoke
end
