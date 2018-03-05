require 'auc_infos'

def species_with_currated_motifs(tf)
  ['mono', 'di'].flat_map{|model_type|
    species_with_currated_motifs_for_model_type(tf, model_type)
  }.uniq
end

def species_with_currated_motifs_for_model_type(tf, model_type)
  currated_motifs = Dir.glob("curation/slices4bench_#{model_type}/#{tf}.*.txt").flat_map{|same_tf_slice_fn|
    File.readlines(same_tf_slice_fn).map(&:strip).reject(&:empty?)
  }.uniq.sort
  species_with_currated_motifs = currated_motifs.select{|motif|
    motif.split('.')[1].start_with?('PEAKS') # ChIP-seq, not legacy hocomoco or other sources
  }.map{|motif|
    motif.split('.').first.split('_').last
  }.uniq
end

['mono', 'di'].each do |model_type|
  task "wlogaucs_for_slices_#{model_type}" do
    all_aucs_all_species = AUCs.all_logaucs_in_folder("auc/#{model_type}/*_datasets/*")

    all_aucs_by_species = {
      'HUMAN' => AUCs.all_logaucs_in_folder("auc/#{model_type}/HUMAN_datasets/*"),
      'MOUSE' => AUCs.all_logaucs_in_folder("auc/#{model_type}/MOUSE_datasets/*"),
    }

    ['HUMAN', 'MOUSE'].each do |target_species|
      FileUtils.mkdir_p "wlogauc/#{model_type}/#{target_species}"
      other_species = ['HUMAN', 'MOUSE'].detect{|s| s != target_species }
      all_aucs_target_species = all_aucs_by_species[target_species]
      all_aucs_other_species = all_aucs_by_species[other_species]
      Dir.glob("curation/slices4bench_#{model_type}/*.txt").each{|fn|
        tf = File.basename(fn, '.txt').split('.').first
        species_to_consider = species_with_currated_motifs(tf).sort

        auc_infos_target_species = AUCs.auc_infos_for_slice(all_aucs_target_species, fn, model_type)
        auc_infos_other_species = AUCs.auc_infos_for_slice(all_aucs_other_species, fn, model_type)
        auc_infos_all_species = AUCs.auc_infos_for_slice(all_aucs_all_species, fn, model_type)
        slice_type = File.basename(fn, '.txt').split('.').last[0]

        consider_target = !auc_infos_target_species.datasets.empty? && species_to_consider.include?(target_species)
        consider_other = !auc_infos_other_species.datasets.empty? && species_to_consider.include?(other_species)


        next  if !consider_target && !consider_other

	# di-processing doesn't know that mono models of target species exist
        has_target_species_models = auc_infos_all_species.models.any?{|model| model.split('.')[0].split('_').last == target_species }
        next  if auc_infos_target_species.datasets.empty? && !has_target_species_models

        if slice_type == 'T'
          semiuniprot = File.basename(fn, '.txt').split('.').first
          uniprot = semiuniprot + '_' + target_species
          next if !auc_infos_target_species.has_only_hocomoco_models? && !auc_infos_target_species.has_chipseq_models?(uniprot, model_type)
        end

        if consider_target && consider_other
          infos = auc_infos_all_species.models.map{|model|
            wauc_target_species = auc_infos_target_species.weighted_auc(model){|dataset| auc_infos_all_species.dataset_quality(dataset) }
            wauc_other_species = auc_infos_other_species.weighted_auc(model){|dataset| auc_infos_all_species.dataset_quality(dataset) }

            wauc = (wauc_target_species * auc_infos_target_species.datasets.size + wauc_other_species).to_f / (auc_infos_target_species.datasets.size + 1)

            best_auc_target_species = auc_infos_target_species.best_auc(model)
            best_auc_other_species = auc_infos_other_species.best_auc(model)
            best_auc = [best_auc_target_species, best_auc_other_species].max
            [model.full_name, wauc, best_auc]
          }
        elsif consider_target
          infos = auc_infos_all_species.models.map{|model|
            wauc = auc_infos_target_species.weighted_auc(model){|dataset| auc_infos_all_species.dataset_quality(dataset) }
            best_auc = auc_infos_target_species.best_auc(model)
            [model.full_name, wauc, best_auc]
          }
        elsif consider_other
          infos = auc_infos_all_species.models.map{|model|
            wauc = auc_infos_target_species.weighted_auc(model){|dataset| auc_infos_all_species.dataset_quality(dataset) }
            best_auc = auc_infos_other_species.best_auc(model)
            [model.full_name, wauc, best_auc]
          }
        else
          $stderr.puts "Check file #{fn} for species #{target_species}"
          infos = []
        end

        infos = infos.sort_by{|model, wauc, maxauc|
          wauc
        }.reverse
        File.write("wlogauc/#{model_type}/#{target_species}/#{File.basename(fn)}", infos.map{|l| l.join("\t") }.join("\n")) unless infos.empty?
      }
    end
  end
end

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
