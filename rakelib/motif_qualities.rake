PCM_EXT = {'mono' => 'pcm', 'di' => 'dpcm'}
PWM_EXT = {'mono' => 'pwm', 'di' => 'dpwm'}
COLLECTION_NAME = {'mono' => 'H11MO', 'di' => 'H11DI'}

# P53_HUMAN.PEAKS033173.gem.control --> P53_HUMAN.PEAKS033173
def peakset(dataset)
  dataset[/^(.*\.PEAKS\d+)\./, 1]
end

def good_peaksets(auc_by_ds, &condition)
  raise 'Specify condition'  unless block_given?
  auc_by_ds.select{|ds, auc|
    condition.call(ds, auc)
  }.map{|ds, auc|
    peakset(ds)
  }.uniq
end

def quality_mark(auc_by_ds, original_species)
  optimal_auc = 0.8
  minimal_auc = 0.65
  auc_by_origin_species_ds = auc_by_ds.select{|ds, auc|
    ds_species = ds.split('.').first.split('_').last
    ds_species == original_species
  }

  peaksets_min = good_peaksets(auc_by_ds){|ds, auc| auc >= minimal_auc }
  peaksets_opt = good_peaksets(auc_by_ds){|ds, auc| auc >= optimal_auc }

  self_dataset = ->(ds){ ds.match(/_#{original_species}.PEAKS/) }
  has_self_dataset_opt = peaksets_opt.any?(&self_dataset)
  has_self_dataset_min = peaksets_min.any?(&self_dataset)

  # At least one peak-set should pass threshold on the peak-set of original species
  # But if it does, we count number of peak-sets passed that threshold over all datasets (both HUMAN and MOUSE)
  # If it doesn't, we assign quality one level less
  if peaksets_opt.count >= 2
    has_self_dataset_opt ? 'A' : 'B'
  elsif peaksets_opt.count >= 1  &&  peaksets_min.count >= 2
    has_self_dataset_min ? 'B' : 'C'
  elsif peaksets_opt.count >= 1  ||  peaksets_min.count >= 2
    has_self_dataset_min ? 'C' : 'D'
  elsif peaksets_min.count >= 1
    has_self_dataset_min ? 'D' : 'E'
  else
    'E'
  end
end


def aucs_from_file(filename)
  return []  unless File.exist?(filename)
  File.readlines(filename).map{|l|
    ds, auc, logauc = l.chomp.split("\t")
    [ds, Float(auc)]
  }
end

def aucs_for_model(model, model_kind)
  [
    "auc/#{model_kind}/HUMAN_datasets/#{model}.txt",
    "auc/#{model_kind}/MOUSE_datasets/#{model}.txt"
  ].map{|auc_fn|
    aucs_from_file(auc_fn)
  }.inject([], &:+)
end

def put_motifs_to_final(infos) #model, final_name, model_kind, should_reverse
  if infos[:model_kind] == 'mono'
    pcm = Bioinform::MotifModel::PCM.from_file(infos[:original_pcm_fn])
    pwm = Bioinform::MotifModel::PWM.from_file(infos[:original_pwm_fn])
  else
    pcm = ModelKind::Di.new.read_pcm(infos[:original_pcm_fn])
    pwm = ModelKind::Di.new.read_pwm(infos[:original_pwm_fn])
  end
  final_name = final_name(infos)
  File.write("final_collection/#{infos[:model_kind]}/pcm/#{final_name}.#{PCM_EXT[infos[:model_kind]]}", (infos[:should_reverse] ? pcm.revcomp : pcm).to_s)
  File.write("final_collection/#{infos[:model_kind]}/pwm/#{final_name}.#{PWM_EXT[infos[:model_kind]]}", (infos[:should_reverse] ? pwm.revcomp : pwm).to_s)
end

def best_model_in_slice(slice_fn)
  model, logauc = File.readlines(slice_fn).map{|l|
    model, logauc, maxlogauc = l.chomp.split("\t")
    [model, Float(logauc)]
  }.max_by{|model, logauc|
    logauc
  }
  # # foolproof check
  # model_second_way = File.readlines(slice_fn).first.split("\t").first
  # raise 'Inconsistent best model'  unless model == model_second_way
  [model, logauc]
end

def collect_novel_motifs(model_kind, species)
  to_reverse = File.readlines('curation/revme_fin.txt').map(&:chomp)
  Dir.glob("wauc/#{model_kind}/#{species}/*.txt").sort.map{|slice_fn|
    model, logauc = best_model_in_slice(slice_fn)
    auc_by_ds = aucs_for_model(model, model_kind)
    slice = File.basename(slice_fn, '.txt')
    quality = quality_mark(auc_by_ds, species)
    [slice, model, logauc, quality]
  }.group_by{|slice, model, logauc, quality|
    slice.split('.').first # semiuniprot
  }.flat_map{|semiuniprot, slices|
    slices.sort_by{|slice, model, logauc, quality|
      # slice_type = slice.split('.').last[0]
      # ['M','S','T'].index(slice_type) # ToDo: сделать порядок по logWAUC
      -logauc
    }.each_with_index.map{|(slice, original_motif, logauc, quality), motif_index|
      uniprot = "#{semiuniprot}_#{species}"
      #
      {
        original_motif: original_motif,
        model_kind: model_kind, species: species,
        uniprot: uniprot, quality: quality, motif_index: motif_index,
        novelty: 'novel', logauc: logauc,
        should_reverse: to_reverse.include?(original_motif.split('~').last),
        original_pcm_fn: "models/pcm/#{model_kind}/all/#{original_motif.split('~')[0]}/#{original_motif}.#{PCM_EXT[model_kind]}",
        original_pwm_fn: "models/pwm/#{model_kind}/all/#{original_motif.split('~')[0]}/#{original_motif}.#{PWM_EXT[model_kind]}",
      }
    }
  }
end

def collect_inherited_motif_infos(inherited_motifs, model_kind)
  inherited_motifs.group_by{|original_motif|
    original_motif.split('.').first # group by uniprot
  }.flat_map{|uniprot, original_motifs|
    species = uniprot.split('_').last
    main_model_quality = original_motifs.map{|original_motif| original_motif.split('.').last }.reject{|quality| quality == 'S' }.first

    original_motifs.sort_by{|original_motif|
      original_motif.split('.').last # quality
    }.each_with_index.map{|original_motif, motif_index|
      # inherited S-models take quality one grade less than the main motif
      original_quality = original_motif.split('.').last
      novel_quality = (original_quality != 'S') ? original_quality : (main_model_quality.ord + 1).chr
      novel_quality = 'D'  if uniprot.start_with?('GLI2_') # manual curation

      {
        original_motif: original_motif,
        model_kind: model_kind, species: species,
        uniprot: uniprot, quality: novel_quality, motif_index: motif_index,
        novelty: 'inherited', logauc: 0,
        should_reverse: false,
        original_pcm_fn: "hocomoco10/#{species}/#{model_kind}/pcm/#{original_motif}.#{PCM_EXT[model_kind]}",
        original_pwm_fn: "hocomoco10/#{species}/#{model_kind}/pwm/#{original_motif}.#{PWM_EXT[model_kind]}",
      }
    }
  }
end

def final_name(infos)
  "#{infos[:uniprot]}.#{COLLECTION_NAME[infos[:model_kind]]}.#{infos[:motif_index]}.#{infos[:quality]}"
end

task 'print_motif_qualities' do
  ['mono', 'di'].flat_map do |model_kind|
    FileUtils.mkdir_p "final_collection/#{model_kind}/pcm/"
    FileUtils.mkdir_p "final_collection/#{model_kind}/pwm/"
    FileUtils.mkdir_p "final_collection/#{model_kind}/logo/"
    novel_motifs = ['HUMAN', 'MOUSE'].flat_map{|species|
      collect_novel_motifs(model_kind, species)
    }

    novel_uniprots = novel_motifs.map{|infos| final_name(infos) }.map{|final_name| final_name.split('.').first }.uniq
    hocomoco10_motifs = Dir.glob("hocomoco10/*/#{model_kind}/pcm/*").map{|pcm_fn|
      File.basename(pcm_fn, File.extname(pcm_fn))
    }

    inherited_motifs = hocomoco10_motifs.reject{|original_motif|
      uniprot = original_motif.split('.').first
      novel_uniprots.include?(uniprot)
    }

    hocomoco_10_infos = collect_inherited_motif_infos(inherited_motifs, model_kind)

    cross_species_infos = []

    infos = novel_motifs + hocomoco_10_infos + cross_species_infos

    motifs_to_ban = ['ERF', 'ETV2_HUMAN', 'MNT_HUMAN\.H10MO\.D', 'MUSC_HUMAN\.H10MO\.D', 'SMRC1', 'ZNF639', 'CLOCK_.*\.H10MO', 'PKNX2', 'YBX1', 'KAISO_MOUSE\.H10MO\.B']
    infos.reject!{|info| # inherit, final_name, model, img
      motifs_to_ban.any?{|motif_pattern|
        pattern = /^(#{motif_pattern}\b|#{motif_pattern}_)/
        info[:uniprot].match(pattern) || info[:original_motif].match(pattern)
      }
    }

    infos.each{|info|
      put_motifs_to_final(info)
    }

    table = infos.map{|infos|
      original_motif, novelty, model_kind, motif_index = infos.values_at(:original_motif, :novelty, :model_kind, :motif_index)
      final_name = final_name(infos)
      [novelty, final_name, original_motif, "<img src='#{model_kind}/logo/#{final_name}_direct.png'>",]
    }

    table.sort_by{|inherit, final_name, model, img|
      final_name
    }.chunk(&:itself).each_slice(300).map{|slice| slice.flat_map(&:last) }
    .each_with_index do |slice, slice_index|
      File.open("final_collection/#{model_kind}_slice_#{slice_index}.html", 'w'){|fw|
        fw.puts "<html><head><style>img{ height:50px; }\ntable,tr,td{ border:1px solid black; }\ntd:first-child{font-weight:bolder;}</style></head><body><table>"
        slice.each do |row|
          fw.puts("<tr>" + row.map{|cell| "<td>#{cell}</td>" }.join + "</tr>")
        end
        fw.puts "</table></body></html>"
      }
    end
  end
end
