require 'auc_infos'
require 'motif_slice'
require 'set'

def species_with_currated_motifs(semiuniprot)
  ['mono', 'di'].flat_map{|model_type|
    species_with_currated_motifs_for_model_type(semiuniprot, model_type)
  }.uniq
end

def species_with_currated_motifs_for_model_type(semiuniprot, model_type)
  Dir.glob("curation/slices4bench_#{model_type}/#{semiuniprot}.*.txt").flat_map{|same_tf_slice_fn|
    MotifsSlice.from_file(same_tf_slice_fn, model_type).species_with_currated_motifs
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


def process_single_slice(motifs_slice, target_species, other_species, target_species_aucs, other_species_aucs, all_species_aucs, model_type)
  species_to_consider = species_with_currated_motifs(motifs_slice.semiuniprot).sort

  dataset_weighting = ->(dataset){ all_species_aucs.dataset_quality(dataset) }

  consider_target = !target_species_aucs.datasets.empty? && species_to_consider.include?(target_species)
  consider_other = !other_species_aucs.datasets.empty? && species_to_consider.include?(other_species)

  return  if !consider_target && !consider_other

  # ToDo: di-processing doesn't know that mono models of target species exist
  has_target_species_models = all_species_aucs.models.any?{|model| model.species == target_species }
  return  if target_species_aucs.datasets.empty? && !has_target_species_models

  if motifs_slice.slice_type == 'T'
    semiuniprot = motifs_slice.semiuniprot
    uniprot = semiuniprot + '_' + target_species
    return  if !target_species_aucs.has_only_hocomoco_models? && !target_species_aucs.has_chipseq_models?(uniprot, model_type)
  end

  if consider_target && consider_other
    infos = all_species_aucs.models.map{|model|
      wauc_target_species = target_species_aucs.weighted_auc(model, &dataset_weighting)
      wauc_other_species = other_species_aucs.weighted_auc(model, &dataset_weighting)

      wauc = (wauc_target_species * target_species_aucs.datasets.size + wauc_other_species).to_f / (target_species_aucs.datasets.size + 1)

      best_auc_target_species = target_species_aucs.best_auc(model)
      best_auc_other_species = other_species_aucs.best_auc(model)
      best_auc = [best_auc_target_species, best_auc_other_species].max
      [model.full_name, wauc, best_auc]
    }
  elsif consider_target
    infos = all_species_aucs.models.map{|model|
      wauc = target_species_aucs.weighted_auc(model, &dataset_weighting)
      best_auc = target_species_aucs.best_auc(model)
      [model.full_name, wauc, best_auc]
    }
  elsif consider_other
    infos = all_species_aucs.models.map{|model|
      wauc = target_species_aucs.weighted_auc(model, &dataset_weighting)
      best_auc = other_species_aucs.best_auc(model)
      [model.full_name, wauc, best_auc]
    }
  else
    $stderr.puts "Check file #{fn} for species #{target_species}"
    infos = []
  end
  return infos
end


['mono', 'di'].each do |model_type|
  task "wlogaucs_for_slices_#{model_type}" do
    semiuniprots = Dir.glob("auc/#{model_type}/*_datasets/*").map{|fn|
      File.basename(fn)
    }.map{|basename|
      basename.split('~').first.split('_').first
    }.uniq.sort

    semiuniprots.each do |semiuniprot|
      all_species_all_aucs = AUCs.in_folder("auc/#{model_type}/*_datasets/#{semiuniprot}_{HUMAN,MOUSE}~*")[:logauc]
      all_aucs_by_species = {
        'HUMAN' => all_species_all_aucs.filter_datasets_by_species('HUMAN'),
        'MOUSE' => all_species_all_aucs.filter_datasets_by_species('MOUSE'),
      }
      ['HUMAN', 'MOUSE'].each do |target_species|
        FileUtils.mkdir_p "wlogauc/#{model_type}/#{target_species}"
        other_species = ['HUMAN', 'MOUSE'].detect{|s| s != target_species }
        Dir.glob("curation/slices4bench_#{model_type}/#{semiuniprot}.*.txt").sort.each{|fn|
          motifs_slice = MotifsSlice.from_file(fn, model_type)

          species_to_consider = species_with_currated_motifs(motifs_slice.semiuniprot).sort

          infos = process_single_slice(
            motifs_slice,
            target_species, other_species,
            all_aucs_by_species[target_species].slice_by_motifs(motifs_slice),
            all_aucs_by_species[other_species].slice_by_motifs(motifs_slice),
            all_species_all_aucs.slice_by_motifs(motifs_slice),
            model_type
          )
          next if !infos || infos.empty?
          infos = infos.sort_by{|model, wauc, maxauc|
            wauc
          }.reverse

          File.write("wlogauc/#{model_type}/#{target_species}/#{File.basename(fn)}", infos.map{|l| l.join("\t") }.join("\n"))  if infos && !infos.empty?
        }
      end

    end
  end
end
