require 'auc_infos'

task :waucs_for_slices do
  ['HUMAN', 'MOUSE'].each do |species|
    FileUtils.mkdir_p "wauc/mono/#{species}"
    all_aucs = AUCs.all_aucs_in_folder("auc/mono/#{species}_datasets/*")
    Dir.glob('curation/slices4bench/*.txt').each{|fn|
      # next  unless File.readlines(fn).any?{|motif_name| motif_name.match(/_#{species}\./) && ! motif_name.match(/\.H10(MO|DI)\./) }
      auc_infos = AUCs.auc_infos_for_slice(all_aucs, fn)
                       .filter_datasets_by_species(species)
      infos = auc_infos.models.map{|model|
        wauc = auc_infos.weighted_auc(model)
        best_auc = auc_infos.best_auc(model)
        [model.full_name, wauc, best_auc ]
      }.sort_by{|model, wauc, maxauc|
        wauc
      }.reverse
      File.write("wauc/mono/#{species}/#{File.basename(fn)}", infos.map{|l| l.join("\t") }.join("\n"))
    }
  end
end

task :dataset_waucs_for_slices do
  ['HUMAN', 'MOUSE'].each do |species|
    FileUtils.mkdir_p "wauc_datasets/#{species}"
    all_aucs = AUCs.all_aucs_in_folder("auc/mono/#{species}_datasets/*")
    Dir.glob('curation/slices4bench/*.txt').each{|fn|
      auc_infos = AUCs.auc_infos_for_slice(all_aucs, fn)
                       .filter_datasets_by_species(species)
      infos = auc_infos.datasets.map{|dataset|
        wauc_ds = auc_infos.dataset_quality(dataset)
        [dataset, wauc_ds]
      }.sort_by{|ds, wauc_ds|
        wauc_ds
      }.reverse
      File.write("wauc_datasets/#{species}/#{File.basename(fn)}", infos.map{|l| l.join("\t") }.join("\n"))
    }
  end
end

