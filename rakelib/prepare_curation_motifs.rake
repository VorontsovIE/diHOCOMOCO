require 'auc_infos'

['mono', 'di'].each do |model_type|
  task "waucs_for_slices_#{model_type}" do
    all_aucs_all_species = AUCs.all_aucs_in_folder("auc/#{model_type}/*_datasets/*")

    all_aucs_by_species = {
      'HUMAN' => AUCs.all_aucs_in_folder("auc/#{model_type}/HUMAN_datasets/*"),
      'MOUSE' => AUCs.all_aucs_in_folder("auc/#{model_type}/MOUSE_datasets/*"),
    }

    ['HUMAN', 'MOUSE'].each do |species|
      FileUtils.mkdir_p "wauc/#{model_type}/#{species}"
      all_aucs = all_aucs_by_species[species]
      other_species = ['HUMAN', 'MOUSE'].detect{|s| s != species }
      all_aucs_other_species = all_aucs_by_species[other_species]
      Dir.glob("curation/slices4bench_#{model_type}/*.txt").each{|fn|
        auc_infos = AUCs.auc_infos_for_slice(all_aucs, fn, model_type)
        auc_infos_other_species = AUCs.auc_infos_for_slice(all_aucs_other_species, fn, model_type)
        auc_infos_all_species = AUCs.auc_infos_for_slice(all_aucs_all_species, fn, model_type)
        slice_type =  File.basename(fn, '.txt').split('.').last[0]

        if slice_type == 'T'
          req_uniprot =  File.basename(fn, '.txt').split('.').first + '_' + species
          chipseq_dataset_code = (model_type == 'mono') ? 'CM' : 'CD'
          if !auc_infos.models.all?{|model| model.full_name.match(/~(DI)?HL~/) }
            if auc_infos.models.none?{|model| model.full_name.start_with?("#{req_uniprot}~#{chipseq_dataset_code}~") }
              next
            end
          end

        end

        infos = auc_infos.models.map{|model|
          if !auc_infos_other_species.datasets.empty?
            wauc_on_species = auc_infos.weighted_auc(model){|dataset| auc_infos_all_species.dataset_quality(dataset) }
            wauc_other_species = auc_infos_other_species.weighted_auc(model){|dataset| auc_infos_all_species.dataset_quality(dataset) }

            wauc = (wauc_on_species * auc_infos.datasets.size + wauc_other_species).to_f / (auc_infos.datasets.size + 1)

            best_auc_on_species = auc_infos.best_auc(model)
            best_auc_other_species = auc_infos_other_species.best_auc(model)
            best_auc = [best_auc_on_species, best_auc_other_species].max
          else
            wauc = auc_infos.weighted_auc(model){|dataset| auc_infos_all_species.dataset_quality(dataset) }
            best_auc_on = auc_infos.best_auc(model)
          end

          [model.full_name, wauc, best_auc ]
        }.sort_by{|model, wauc, maxauc|
          wauc
        }.reverse
        File.write("wauc/#{model_type}/#{species}/#{File.basename(fn)}", infos.map{|l| l.join("\t") }.join("\n")) unless infos.empty?
      }
    end
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

