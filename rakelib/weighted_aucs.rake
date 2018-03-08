require 'auc_infos'
require 'set'

def species_with_currated_motifs(tf)
  ['mono', 'di'].flat_map{|model_type|
    species_with_currated_motifs_for_model_type(tf, model_type)
  }.uniq
end

def species_with_currated_motifs_for_model_type(tf, model_type)
  Dir.glob("curation/slices4bench_#{model_type}/#{tf}.*.txt").flat_map{|same_tf_slice_fn|
    MotifSlice.from_file(same_tf_slice_fn, model_type).species_with_currated_motifs
  }.uniq.sort
end

def tf_for_species_exist?(semiuniprot, species)
  $species_for_tf ||= begin
    result = Hash.new{|h,k| h[k] = Set.new }
    Dir.glob('models/pwm/*/all/*/*.{d,}pwm').map{|fn|
      File.basename(fn, '.pwm')
    }.each{|motif|
      semiuniprot, species = motif.split('~').first.split('_')
      result[semiuniprot] << species
    }

    Dir.glob('control/control/*.mfa').map{|fn|
      File.basename(fn, '.mfa')
    }.each{|dataset|
      semiuniprot, species = dataset.split('.').first.split('_')
      result[semiuniprot] << species
    }
    result
  end
  raise "Unknown TF `#{semiuniprot}` for any species"  if $species_for_tf[semiuniprot].empty?
  $species_for_tf[semiuniprot].include?(species)
end

['mono', 'di'].each do |model_type|
  task "wlogaucs_for_slices_#{model_type}" do
    all_aucs_all_species = AUCs.in_folder("auc/#{model_type}/*_datasets/*")[:logauc]

    all_aucs_by_species = {
      'HUMAN' => AUCs.in_folder("auc/#{model_type}/HUMAN_datasets/*")[:logauc],
      'MOUSE' => AUCs.in_folder("auc/#{model_type}/MOUSE_datasets/*")[:logauc],
    }

    ['HUMAN', 'MOUSE'].each do |target_species|
      FileUtils.mkdir_p "wlogauc/#{model_type}/#{target_species}"
      other_species = ['HUMAN', 'MOUSE'].detect{|s| s != target_species }
      all_aucs_target_species = all_aucs_by_species[target_species]
      all_aucs_other_species = all_aucs_by_species[other_species]
      Dir.glob("curation/slices4bench_#{model_type}/*.txt").each{|fn|

        motifs_slice = MotifsSlice.from_file(fn, model_type)

        tf = motifs_slice.semiuniprot
        slice_type = motifs_slice.slice_type
        species_to_consider = species_with_currated_motifs(tf).sort

        auc_infos_target_species = all_aucs_target_species.slice(motifs_slice)
        auc_infos_other_species = all_aucs_other_species.slice(motifs_slice)
        auc_infos_all_species = all_aucs_all_species.slice(motifs_slice)
        dataset_weighting = ->(dataset){ auc_infos_all_species.dataset_quality(dataset) }

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
            wauc_target_species = auc_infos_target_species.weighted_auc(model, &dataset_weighting)
            wauc_other_species = auc_infos_other_species.weighted_auc(model, &dataset_weighting)

            wauc = (wauc_target_species * auc_infos_target_species.datasets.size + wauc_other_species).to_f / (auc_infos_target_species.datasets.size + 1)

            best_auc_target_species = auc_infos_target_species.best_auc(model)
            best_auc_other_species = auc_infos_other_species.best_auc(model)
            best_auc = [best_auc_target_species, best_auc_other_species].max
            [model.full_name, wauc, best_auc]
          }
        elsif consider_target
          infos = auc_infos_all_species.models.map{|model|
            wauc = auc_infos_target_species.weighted_auc(model, &dataset_weighting)
            best_auc = auc_infos_target_species.best_auc(model)
            [model.full_name, wauc, best_auc]
          }
        elsif consider_other
          infos = auc_infos_all_species.models.map{|model|
            wauc = auc_infos_target_species.weighted_auc(model, &dataset_weighting)
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
