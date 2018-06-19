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
  Rake::Task['collect_and_normalize_data'].invoke

  Rake::Task['calculate_local_backgrounds'].invoke
  Rake::Task['average_local_backgrounds'].invoke

  Rake::Task['precalculate_thresholds'].invoke
  Rake::Task['group_controls'].invoke
  Rake::Task['calculate_occurence_scores'].invoke

  Rake::Task['choose_motifs_for_final_collection'].invoke
  Rake::Task['put_motifs_into_final_collection'].invoke
  Rake::Task['sequence_logos'].invoke
  Rake::Task['precalc_thresholds_for_final_bundle'].invoke
  # Rake::Task['put_thresholds_to_json'].invoke # optional
  Rake::Task['repack_final_collection'].invoke
  ## --> Don't forget to manually put comments to annotations here
  Rake::Task['family_tree_svg'].invoke
  Rake::Task['archive_final_bundle'].invoke
end
