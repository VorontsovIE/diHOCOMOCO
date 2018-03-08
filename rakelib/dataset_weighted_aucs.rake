require 'auc_infos'

task :dataset_waucs_for_slices do
 ['mono', 'di'].each do |model_type|
    FileUtils.mkdir_p "wlogauc_datasets/#{model_type}"
    all_aucs = auc_infos_in_folder("auc/#{model_type}/*_datasets/*")
    Dir.glob("curation/slices4bench_#{model_type}/*.txt").each{|fn|
      motifs_slice = MotifsSlice.from_file(fn, model_type)
      auc_infos = all_aucs[:auc].auc_infos_for_slice(motifs_slice)
      logauc_infos = all_aucs[:logauc].auc_infos_for_slice(motifs_slice)
      datasets = auc_infos.datasets
      infos = datasets.map{|ds|
        [ds, auc_infos.dataset_quality(ds), logauc_infos.dataset_quality(ds)]
      }.sort_by{|ds, ds_wauc, ds_wlogauc| ds_wlogauc }.reverse
      File.write("wauc_datasets/#{model_type}/#{File.basename(fn)}", infos.map{|l| l.join("\t") }.join("\n"))
    }
  end
end
