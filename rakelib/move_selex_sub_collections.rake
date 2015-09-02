def move_to_folder(glob, dest_folder)
  FileList[glob].each{|fn|
    mkdir_p  File.join(dest_folder, File.dirname(fn))
    mv fn, File.join(dest_folder, fn)
  }
end

desc 'Move collections SMS/SDS to a separate folder: we don\'t use them anymore. Task made not to recalculate results of computationally expensive steps'
task :move_selex_sub_collections do
  mkdir_p 'removed_selex_sub'

  move_to_folder("models/pcm/mono/all/*/*~SMS~*.pcm", 'removed_selex_sub')
  move_to_folder("models/pwm/mono/all/*/*~SMS~*.pwm", 'removed_selex_sub')
  move_to_folder("models/thresholds/mono/all/*/*~SMS~*.thr", 'removed_selex_sub')

  move_to_folder("models/pcm/di/all/*/*~SDS~*.dpcm", 'removed_selex_sub')
  move_to_folder("models/pwm/di/all/*/*~SDS~*.dpwm", 'removed_selex_sub')
  move_to_folder("models/thresholds/di/all/*/*~SDS~*.thr", 'removed_selex_sub')

  ['SMS','SDS'].each do |collection|
    move_to_folder("models/logo/*/*~#{collection}~*.png", 'removed_selex_sub')
  end

  move_to_folder("occurences/scores/mono/*/*~SMS~*/", 'removed_selex_sub')
  move_to_folder("occurences/scores/di/*/*~SDS~*/", 'removed_selex_sub')

  move_to_folder("occurences/pvalues/*/*~SMS~*/", 'removed_selex_sub')
  move_to_folder("occurences/pvalues/*/*~SDS~*/", 'removed_selex_sub')

  move_to_folder("occurences/corrected_pvalues/*/*~SMS~*/", 'removed_selex_sub')
  move_to_folder("occurences/corrected_pvalues/*/*~SDS~*/", 'removed_selex_sub')

  move_to_folder("occurences/auc/*/*~SMS~*.txt", 'removed_selex_sub')
  move_to_folder("occurences/auc/*/*~SDS~*.txt", 'removed_selex_sub')
end
