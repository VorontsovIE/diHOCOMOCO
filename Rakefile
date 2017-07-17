$:.unshift File.absolute_path('lib')
require 'rake/clean'

desc 'Move all motifs/words into their folders, normalize their names, convert to PWMs and so on'
task 'collect_and_normalize_data' do
  Rake::Task['collect_and_normalize_data:collect_pcm'].invoke
  Rake::Task['collect_and_normalize_data:rename_motifs'].invoke
  Rake::Task['collect_and_normalize_data:convert_pcm_to_pwm'].invoke
  Rake::Task['collect_and_normalize_data:collect_words'].invoke
  Rake::Task['collect_and_normalize_data:rename_words'].invoke
end

task :default do
  Rake::Task['unpack'].invoke # 4m
  Rake::Task['collect_and_normalize_data'].invoke # 1m

  Rake::Task['calculate_local_backgrounds'].invoke # 5m
  Rake::Task['average_local_backgrounds'].invoke

  Rake::Task['precalculate_thresholds'].invoke # 8950m; parallelizable

  Rake::Task['calculate_occurence_scores'].invoke

  Rake::Task['sequence_logos'].invoke # 30m; parallelizable
  Rake::Task['make_final_collection'].invoke
  Rake::Task['final_collection_summary'].invoke
end
