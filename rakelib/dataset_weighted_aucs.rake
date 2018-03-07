require 'auc_infos'

task :dataset_wlogaucs_for_slices do
 ['mono', 'di'].each do |model_type|
    FileUtils.mkdir_p "wlogauc_datasets/#{model_type}"
    all_aucs = AUCs.all_logaucs_in_folder("auc/#{model_type}/*_datasets/*")
    Dir.glob("curation/slices4bench_#{model_type}/*.txt").each{|fn|
      auc_infos = AUCs.auc_infos_for_slice(all_aucs, fn, model_type)
      infos = auc_infos.datasets.map{|dataset|
        wauc_ds = auc_infos.dataset_quality(dataset)
        [dataset, wauc_ds]
      }.sort_by{|ds, wauc_ds|
        wauc_ds
      }.reverse
      File.write("wlogauc_datasets/#{model_type}/#{File.basename(fn)}", infos.map{|l| l.join("\t") }.join("\n"))
    }
  end
end

task :dataset_waucs_for_slices do
 ['mono', 'di'].each do |model_type|
    FileUtils.mkdir_p "wauc_datasets/#{model_type}"
    all_aucs = AUCs.all_aucs_in_folder("auc/#{model_type}/*_datasets/*")
    Dir.glob("curation/slices4bench_#{model_type}/*.txt").each{|fn|
      auc_infos = AUCs.auc_infos_for_slice(all_aucs, fn, model_type)
      infos = auc_infos.datasets.map{|dataset|
        wauc_ds = auc_infos.dataset_quality(dataset)
        [dataset, wauc_ds]
      }.sort_by{|ds, wauc_ds|
        wauc_ds
      }.reverse
      File.write("wauc_datasets/#{model_type}/#{File.basename(fn)}", infos.map{|l| l.join("\t") }.join("\n"))
    }
  end
end
