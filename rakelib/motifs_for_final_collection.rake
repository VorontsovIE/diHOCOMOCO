require 'set'
require 'json'
require 'uniprot_info'
require 'motif_family_recognizer'

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

# we should get data for both species in order to know origin of motif obtained as cross-species
def original_motifs_origin
  @original_motifs_origin_cache ||= begin
    ['mono', 'di'].flat_map{|arity|
      ['HUMAN', 'MOUSE'].flat_map{|species|
        File.readlines("hocomoco10/#{species}/#{arity}/final_collection.tsv").drop(1).map{|line|
          motif, hocomoco10_motif_origin = line.chomp.split("\t").values_at(0, 12)
          [motif, origin_by_motif_in_hocomoco10(hocomoco10_motif_origin)]
        }
      }
    }.to_h
  end
end

def get_release_and_source(original_motif)
  if original_motif.match(/~(CM|CD)~/)
    ['HOCOMOCOv11', 'ChIP-Seq']
  else
    prev_motif = original_motif.split('~').last
    if original_motifs_origin[prev_motif] == 'HOCOMOCO v9'
      ['HOCOMOCOv9', 'Integrative']
    else
      ['HOCOMOCOv10', original_motifs_origin[prev_motif]]
    end
  end
end

def infos_by_uniprot_id
  @infos_by_uniprot_id_cache ||= begin
    UniprotInfo.each_in_file('uniprot_HomoSapiens_and_MusMusculus_lots_of_infos.tsv')
               .group_by(&:uniprot_id)
  end
end

def motif_words_filename(model_kind, original_motif_name)
  original_name = original_motif_name.split('~').last  # name w/o 'uniprot~collection~' part
  if original_name.match(/\.H10(MO|DI)\./)
    "models/words/#{model_kind}/hocomoco_legacy/#{original_name}.words"
  else
    "models/words/#{model_kind}/chipseq/#{original_name}.#{model_kind.wordlist_extension}"
  end
end

def motif_words(words_fn)
  if File.extname(words_fn) == '.words'
    File.readlines(words_fn).map(&:strip).reject(&:empty?)
  else # .list /.dlist
    File.readlines(words_fn).drop(1).map{|l|
      l.chomp.split("\t")[2]
    }
  end
end

def motif_infos_dump(infos, arity)
  model_kind = ModelKind.get(infos[:model_kind])

  bundle_list = ['full']
  bundle_list << 'core'  if infos[:motif_index] == 0 && ['A','B','C'].include?(infos[:quality])

  pcm = model_kind.read_pcm(infos[:original_pcm_fn])
  pwm = model_kind.read_pwm(infos[:original_pwm_fn])

  words_fn = motif_words_filename(model_kind, infos[:original_motif])
  release, source = get_release_and_source(infos[:original_motif])

  best_auc_human, num_datasets_human = get_auc_stats("auc/#{arity}/HUMAN_datasets/#{infos[:original_motif]}.txt").values_at(:best_auc, :num_datasets)
  best_auc_mouse, num_datasets_mouse = get_auc_stats("auc/#{arity}/MOUSE_datasets/#{infos[:original_motif]}.txt").values_at(:best_auc, :num_datasets)

  uniprot = infos[:uniprot]
  uniprot_infos = infos_by_uniprot_id[uniprot]

  recognizers_by_level = PROTEIN_FAMILY_RECOGNIZERS[infos[:species]]
  motif_superclass = recognizers_by_level[1].subfamilies_by_uniprot_id(uniprot)
  motif_class = recognizers_by_level[2].subfamilies_by_uniprot_id(uniprot)
  motif_families = recognizers_by_level[3].subfamilies_by_uniprot_id(uniprot)
  motif_subfamilies = recognizers_by_level[4].subfamilies_by_uniprot_id(uniprot)
  motif_genus = recognizers_by_level[5].subfamilies_by_uniprot_id(uniprot)

  {
    name: final_name(infos),
    bundle_list: bundle_list,
    should_reverse: infos[:should_reverse],
    comments: "",
    model_kind: infos[:model_kind], # aka `arity`
    original_motif: infos[:original_motif],
    release: release,
    source: source,
    species: infos[:species],
    uniprot: infos[:uniprot],
    quality: infos[:quality],
    motif_index: infos[:motif_index], # aka `rank`
    novelty: infos[:novelty],
    logauc: infos[:logauc],
    length: pcm.length,
    best_auc_human: best_auc_human, num_datasets_human: num_datasets_human,
    best_auc_mouse: best_auc_mouse, num_datasets_mouse: num_datasets_mouse,
    uniprot_infos: {
      uniprot_acs: uniprot_infos.flat_map(&:uniprot_ac),
      gene_names: uniprot_infos.flat_map(&:primary_gene_name),
      hgnc_ids: uniprot_infos.flat_map(&:hgnc_ids),
      mgi_ids: uniprot_infos.flat_map(&:mgi_ids),
      entrezgene_ids: uniprot_infos.flat_map(&:entrezgene_ids),
    },
    tfclass: {
      motif_superclasses: motif_superclass.map(&:to_s),
      motif_classes: motif_class.map(&:to_s),
      motif_families: motif_families.map(&:to_s),
      motif_subfamilies: motif_subfamilies.map(&:to_s),
      motif_geni: motif_genus.map(&:to_s),
    },
    pcm_fn: infos[:original_pcm_fn],
    pwm_fn: infos[:original_pwm_fn],
    words_fn: words_fn,
  }
end

def process_motif_dump(infos)
  model_kind = ModelKind.get(infos[:model_kind])
  final_name = infos[:name]
  pcm = model_kind.create_pcm(infos[:pcm]).named(final_name)
  pwm = model_kind.create_pwm(infos[:pwm]).named(final_name)
  File.write("final_collection/#{model_kind}/pcm/#{final_name}.#{model_kind.pcm_extension}", pcm.to_s)
  File.write("final_collection/#{model_kind}/pwm/#{final_name}.#{model_kind.pwm_extension}", pwm.to_s)
  File.write("final_collection/#{model_kind}/words/#{final_name}.words", infos[:words].map{|x| "#{x}\n" }.join)
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

def final_name(infos)
  "#{infos[:uniprot]}.#{COLLECTION_NAME[infos[:model_kind]]}.#{infos[:motif_index]}.#{infos[:quality]}"
end

def hocomoco10_motifs(model_kind)
  Dir.glob("hocomoco10/*/#{model_kind}/pcm/*").map{|pcm_fn|
    File.basename(pcm_fn, File.extname(pcm_fn))
  }
end

def collect_novel_motifs(model_kind, species)
  to_reverse_mono = File.readlines('curation/to_reverse_mono.txt').map(&:strip)
  to_reverse_di = File.readlines('curation/to_reverse_di.txt').map(&:strip)
  uniprots_failed_curation = File.readlines('curation/uniprots_failed_curation.txt').map(&:strip)
  Dir.glob("wlogauc/#{model_kind}/#{species}/*.txt").sort.map{|slice_fn|
    model, logauc = best_model_in_slice(slice_fn)
    auc_by_ds = aucs_for_model(model, model_kind)
    slice = File.basename(slice_fn, '.txt')
    quality = quality_mark(auc_by_ds, species)
    [slice, model, logauc, quality]
  }.group_by{|slice, model, logauc, quality|
    slice.split('.').first # semiuniprot
  }.flat_map{|semiuniprot, slices|
    slices.sort_by{|slice, model, logauc, quality|
      -logauc
    }.reject{|slice, model, logauc, quality|
      quality >= 'E'
    }.each_with_index.map{|(slice, original_motif, logauc, quality), motif_index|
      uniprot = "#{semiuniprot}_#{species}"
      if model_kind == 'mono'
        should_reverse = to_reverse_mono.include?(original_motif.split('~').last)
      else
        final_name = "#{uniprot}.#{COLLECTION_NAME[model_kind]}.#{motif_index}.#{quality}"
        should_reverse = to_reverse_di.include?(final_name)
      end
      {
        original_motif: original_motif,
        model_kind: model_kind, species: species,
        uniprot: uniprot, quality: quality, motif_index: motif_index,
        novelty: 'novel', logauc: logauc,
        should_reverse: should_reverse,
        original_pcm_fn: "models/pcm/#{model_kind}/all/#{original_motif.split('~')[0]}/#{original_motif}.#{PCM_EXT[model_kind]}",
        original_pwm_fn: "models/pwm/#{model_kind}/all/#{original_motif.split('~')[0]}/#{original_motif}.#{PWM_EXT[model_kind]}",
      }
    }
  }.reject{|info|
    uniprots_failed_curation.include?(info[:uniprot].split('_').first)
  }
end

def inherited_motifs_infos_for_tf(uniprot, hocomoco10_tf_motifs, model_kind)
  to_reverse_mono = File.readlines('curation/to_reverse_mono.txt').map(&:strip)
  to_reverse_di = File.readlines('curation/to_reverse_di.txt').map(&:strip)
  species = uniprot.split('_').last
  main_model_quality = hocomoco10_tf_motifs.map{|original_motif| original_motif.split('.').last }.reject{|quality| quality == 'S' }.first

  hocomoco10_tf_motifs.sort_by{|original_motif|
    original_motif.split('.').last # quality
  }.each_with_index.map{|original_motif, motif_index|
    original_quality = original_motif.split('.').last
    novel_quality = original_quality
    novel_quality = 'D'  if uniprot.start_with?('GLI2_') # manual curation
    novel_quality = main_model_quality.succ  if original_quality == 'S'
    novel_quality = 'D'  if uniprot == 'ZBTB4_HUMAN' # manual curation (not to drop ZBTB4 of S-quality which is a motif for methylated DNA)
    if model_kind == 'mono'
      should_reverse = to_reverse_mono.include?(original_motif.split('~').last)
    else
      final_name = "#{uniprot}.#{COLLECTION_NAME[model_kind]}.#{motif_index}.#{novel_quality}"
      should_reverse = to_reverse_di.include?(final_name)
    end
    {
      original_motif: original_motif,
      model_kind: model_kind, species: species,
      uniprot: uniprot, quality: novel_quality, motif_index: motif_index,
      novelty: 'inherited', logauc: 0,
      should_reverse: should_reverse,
      original_pcm_fn: "hocomoco10/#{species}/#{model_kind}/pcm/#{original_motif}.#{PCM_EXT[model_kind]}",
      original_pwm_fn: "hocomoco10/#{species}/#{model_kind}/pwm/#{original_motif}.#{PWM_EXT[model_kind]}",
    }
  }.reject{|infos|
    infos[:quality] >= 'E'
  }
end

def collect_inherited_motif_infos(inherited_motifs, model_kind)
  inherited_motifs.group_by{|original_motif|
    original_motif.split('.').first # group by uniprot
  }.flat_map{|uniprot, original_tf_motifs|
    inherited_motifs_infos_for_tf(uniprot, original_tf_motifs, model_kind)
  }
end


def other_species(uniprot)
  uniprot.match(/_HUMAN/) ? uniprot.sub(/_HUMAN/, '_MOUSE') : uniprot.sub(/_MOUSE/, '_HUMAN')
end

def cross_species_infos(chipseq_infos, model_kind)
  chipseq_uniprots = chipseq_infos.map{|info| info[:uniprot] }.to_set
  can_be_cross_assigned = chipseq_uniprots.map{|uniprot| other_species(uniprot) }.to_set - chipseq_uniprots
  hocomoco10_motifs = hocomoco10_motifs(model_kind)
  hocomoco10_uniprots = hocomoco10_motifs.map{|motif| motif.split('.').first }.to_set
  should_be_cross_assigned = can_be_cross_assigned & hocomoco10_uniprots
  # puts ['Uniprot', 'origin', 'chosen motifs', 'alternative was'].join("\t")
  should_be_cross_assigned.flat_map{|uniprot|
    hocomoco10_variants = hocomoco10_motifs.select{|motif| motif.split('.').first == uniprot }
    chipseq_variants_infos = chipseq_infos.select{|info| other_species(info[:uniprot]) == uniprot }
    chipseq_variants = chipseq_variants_infos.map{|info| info[:original_motif] }
    hocomoco10_best_quality = hocomoco10_variants.map{|motif| motif.split('.').last }.min
    chipseq_best_quality = chipseq_variants_infos.map{|info| info[:quality] }.min
    chipseq_motifs = chipseq_variants_infos.select{|info| info[:motif_index] < 2 }.map{|info| info[:original_motif] }
    # chipseq_motifs_number_2 = chipseq_variants_infos.select{|info| info[:motif_index] == 2 }.map{|info| info[:original_motif] }
    # if hocomoco10_best_quality < (chipseq_best_quality.ord + 1).chr
    #   puts [uniprot, 'Hocomoco10', hocomoco10_variants.join("; "), (chipseq_motifs + chipseq_motifs_number_2).join("; ")].join("\t")
    # else
    #   puts [uniprot, "Cross-species", chipseq_motifs.join("; "), (hocomoco10_variants + chipseq_motifs_number_2).join("; ")].join("\t")
    # end
    if hocomoco10_best_quality < (chipseq_best_quality.ord + 1).chr
      inherited_motifs_infos_for_tf(uniprot, hocomoco10_variants, model_kind).map{|info|
        info.merge(novelty: 'inherited (better than cross-species)')
      }
    else
      chipseq_variants_infos.map{|info|
        info.merge({
          uniprot: uniprot,
          species: uniprot.split('_').last,
          quality: (info[:quality].ord + 1).chr,
          novelty: 'cross-species',
        })
      }.reject{|info|
        info[:quality] >= 'E'
      }
    end
  }
end

desc 'Select which motifs and with which names are put into collection (create JSON files for motifs)'
task 'choose_motifs_for_final_collection' do
  ['mono', 'di'].each do |model_kind|
    FileUtils.mkdir_p "final_collection/#{model_kind}/json/"
    novel_chipseq_infos = ['HUMAN', 'MOUSE'].flat_map{|species|
      collect_novel_motifs(model_kind, species)
    }

    cross_species_infos = cross_species_infos(novel_chipseq_infos, model_kind)

    novel_semiuniprots = novel_chipseq_infos.map{|infos|
      final_name(infos)
    }.map{|final_name|
      final_name.split('.').first
    }.map{|uniprot|
      uniprot.split('_').first
    }.uniq

    inherited_motifs = hocomoco10_motifs(model_kind).reject{|original_motif|
      uniprot = original_motif.split('.').first
      semiuniprot = uniprot.split('_').first
      novel_semiuniprots.include?(semiuniprot)
    }

    hocomoco10_infos = collect_inherited_motif_infos(inherited_motifs, model_kind)

    motifs_to_ban = ['ERF', 'ETV2_HUMAN', 'MNT_HUMAN\.H10MO\.D', 'MUSC_HUMAN\.H10MO\.D', 'SMRC1', 'ZN639', 'CLOCK_.*\.H10MO', 'PKNX2', 'YBOX1', 'KAISO_MOUSE\.H10MO\.B', 'GABP1']
    hocomoco10_infos.reject!{|info| # inherit, final_name, model, img
      motifs_to_ban.any?{|motif_pattern|
        pattern = /^(#{motif_pattern}\b|#{motif_pattern}_)/
        info[:uniprot].match(pattern) || info[:original_motif].match(pattern)
      }
    }

    infos = novel_chipseq_infos + hocomoco10_infos + cross_species_infos

    novel_chipseq_uniprots = novel_chipseq_infos.map{|info| info[:uniprot] }.to_set
    hocomoco10_uniprots = hocomoco10_infos.map{|info| info[:uniprot] }.to_set
    cross_species_uniprots = cross_species_infos.map{|info| info[:uniprot] }.to_set
    raise 'Overlap'  if novel_chipseq_uniprots.intersect?(hocomoco10_uniprots) || novel_chipseq_uniprots.intersect?(cross_species_uniprots) || cross_species_uniprots.intersect?(hocomoco10_uniprots)
    unless hocomoco10_motifs(model_kind).map{|motif| motif.split('.').first }.all?{|uniprot| infos.map{|info| info[:uniprot] }.include?(uniprot) }
      $stderr.puts "Some HOCOMOCO v10 motifs not covered by new collection:"
      $stderr.puts (hocomoco10_motifs(model_kind).map{|motif| motif.split('.').first } - infos.map{|info| info[:uniprot] }).uniq
    end

    infos.each{|info|
      infos_dump = motif_infos_dump(info, model_kind)
      json_filename = "final_collection/#{infos_dump[:model_kind]}/json/#{infos_dump[:name]}.json"
      File.write(json_filename, infos_dump.to_json)
    }
  end
end

###################
desc 'Take motif JSON files and put them into collection'
task 'put_motifs_into_final_collection' do
  ['mono', 'di'].each do |model_kind|
    motif_infos = Dir.glob("final_collection/#{model_kind}/json_processed/*.json").map{|json_fn|
      JSON.parse(File.read(json_fn), symbolize_names: true)
    }

    FileUtils.mkdir_p "final_collection/#{model_kind}/pcm/"
    FileUtils.mkdir_p "final_collection/#{model_kind}/pwm/"
    FileUtils.mkdir_p "final_collection/#{model_kind}/words/"
    motif_infos.each{|infos|
      # put pcm/pwm/words into final their folders
      process_motif_dump(infos)
    }

    table = motif_infos.map{|infos|
      original_motif, novelty, model_kind, motif_index, final_name = infos.values_at(:original_motif, :novelty, :model_kind, :motif_index, :name)
      [novelty, final_name, original_motif, "#{model_kind}/logo/#{final_name}.png",]
    }

    File.open("final_collection/#{model_kind}.html", 'w') do |fw|
      fw.puts "<html><head><style>img{ height:50px; }\ntable,tr,td{ border:1px solid black; }\ntd:first-child{font-weight:bolder;}</style></head><body><table>"
      table.sort_by{|inherit, final_name, model, img_src|
        final_name
      }.each{|inherit, final_name, model, img_src|
        fw.puts("<tr>" + [inherit, final_name, model, "<img src='#{img_src}'>"].map{|cell| "<td>#{cell}</td>" }.join + "</tr>")
      }
      fw.puts "</table></body></html>"
    end

    File.open("final_collection/#{model_kind}.tsv", 'w'){|fw|
      fw.puts ['origin', 'final_name', 'original_motif', 'img'].join("\t")
      table.sort_by{|inherit, final_name, model, img|
        final_name
      }.each{|row|
        fw.puts(row.join("\t"))
      }
    }
  end
end
