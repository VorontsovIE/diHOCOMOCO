# rake update_tfs["BMAL1_MOUSE GATA1_MOUSE GRHL1_HUMAN NFIC_HUMAN SP1_HUMAN TAL1_HUMAN TAL1_MOUSE FOXH1_HUMAN RCOR1_HUMAN"]
task :update_tfs, :tfs_to_update do |t, args|
  tfs_to_update = args[:tfs_to_update].split

  tfs_to_update.each do |tf|
    FileList["models/thresholds/mono/all/#{tf}^*"].each{|folder| rm_rf folder }
    FileList["models/thresholds/di/all/#{tf}^*"].each{|folder| rm_rf folder }
    rm_rf "occurences/scores/mono/#{tf}"
    rm_rf "occurences/scores/di/#{tf}"
    rm_rf "occurences/pvalues/#{tf}"
    rm_rf "occurences/corrected_pvalues/#{tf}"
    rm_rf "occurences/auc/#{tf}"
  end
end
