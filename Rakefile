$:.unshift File.absolute_path('lib')
require 'rake/clean'

desc 'Move all motifs/words into their folders, normalize their names, convert to PWMs and so on'
task 'collect_and_normalize_data' do
  Rake::Task['collect_and_normalize_data:collect_pcm'].invoke
  Rake::Task['collect_and_normalize_data:rename_motifs'].invoke
  Rake::Task['collect_and_normalize_data:convert_pcm_to_pwm'].invoke
end

task :default do
  Rake::Task['unpack'].invoke # 4m
  Rake::Task['collect_and_normalize_data'].invoke # 1m
  Rake::Task['remove_non_currated_models'].invoke # 17s

  Rake::Task['calculate_local_backgrounds'].invoke # 5m
  Rake::Task['precalculate_thresholds'].invoke # 8950m; parallelizable

  Rake::Task['calculate_occurence_scores'].invoke # 770m; parallelizable
  Rake::Task['scores_to_pvalues'].invoke # 3570m; parallelizable; Sometimes can hang, kill hanged process remove result and run again
  Rake::Task['correct_pvalues'].invoke # 5m25s
  Rake::Task['calculate_auc'].invoke # 2m

  Rake::Task['sequence_logos'].invoke # 30m; parallelizable
  Rake::Task['make_final_collection'].invoke
  Rake::Task['final_collection_summary'].invoke
end
