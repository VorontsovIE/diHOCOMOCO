require 'models'
require 'json'
require 'information_content'
require 'auc_infos'
require 'motif_family_recognizer'
require 'ape_find_threshold'
require 'formatters'
require 'thresholds_bsearch'

# whole collection in a single file (one for all PCMs, one for all PWMs etc)
def save_collection_in_single_files!(folder, species, arity, requested_pvalues, thresholds_by_model, hocomoco_prefix)
  model_kind = ModelKind.get(arity)
  pcms = Dir.glob("#{folder}/pcm/*").sort.map{|fn| model_kind.read_pcm(fn) }
  pwms = Dir.glob("#{folder}/pwm/*").sort.map{|fn| model_kind.read_pwm(fn) }

  File.write File.join(folder, "#{hocomoco_prefix}pcms_#{species}_#{arity}.txt"), pcms.map(&:to_s).join("\n")
  File.write File.join(folder, "#{hocomoco_prefix}pwms_#{species}_#{arity}.txt"), pwms.map(&:to_s).join("\n")

  if arity == 'mono'
    File.write File.join(folder, "#{hocomoco_prefix}#{species}_mono_meme_format.meme"), in_meme_format(pcms)
    File.write File.join(folder, "#{hocomoco_prefix}#{species}_mono_transfac_format.txt"), in_transfac_format(pcms)
    File.write File.join(folder, "#{hocomoco_prefix}#{species}_mono_jaspar_format.txt"), in_jaspar_format(pcms)
    requested_pvalues.each do |requested_pvalue|
      File.write File.join(folder, "#{hocomoco_prefix}#{species}_mono_homer_format_#{requested_pvalue}.motif"), in_homer_format(pcms, thresholds_by_model, pvalue: requested_pvalue)
    end
  end
end

def origin_by_motif_in_hocomoco10(motif)
  collection_name = motif.split('~')[1]
  {
    'CD' => 'ChIP-Seq', # Actually these models aren't in collection
    'CM' => 'ChIP-Seq',
    'PAPAM' => 'ChIP-Seq',
    'PAPAD' => 'ChIP-Seq',

    'SMF' => 'HT-SELEX',
    'SMI' => 'HT-SELEX',
    'SDF' => 'HT-SELEX', # Actually these models aren't in collection
    'SDI' => 'HT-SELEX', # Actually these models aren't in collection

    'HL' => 'HOCOMOCO v9',
  }[collection_name]
end

def get_auc_stats(auc_fn)
  return({ best_auc: nil, num_datasets: nil })  unless File.exist?(auc_fn)
  aucs = File.readlines(auc_fn).map{|l|
    Float(l.chomp.split("\t")[1])
  }
  { best_auc: aucs.max, num_datasets: aucs.size }
end

def get_release_and_source(original_motif, original_motifs_origin)
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

def store_collection_summary(species, arity, motif_names, num_words_by_motif:, pcm_by_motif:, original_motifs:, original_motifs_origin:, infos_by_uniprot_id:, output_stream:)
  recognizers_by_level = PROTEIN_FAMILY_RECOGNIZERS[species]

  model_infos = motif_names.map{|motif|
    pcm = pcm_by_motif[motif]
    num_words = num_words_by_motif[motif]

    original_motif = original_motifs[motif]

    best_auc_human, num_datasets_human = get_auc_stats("auc/#{arity}/HUMAN_datasets/#{original_motif}.txt").values_at(:best_auc, :num_datasets)
    best_auc_mouse, num_datasets_mouse = get_auc_stats("auc/#{arity}/MOUSE_datasets/#{original_motif}.txt").values_at(:best_auc, :num_datasets)

    release, source = get_release_and_source(original_motif, original_motifs_origin)

    uniprot = motif.split('.').first
    quality = motif.split('.').last
    rank = motif.split('.')[-2].to_i
    uniprot_infos = infos_by_uniprot_id[uniprot]
    motif_families = recognizers_by_level[3].subfamilies_by_uniprot_id(uniprot)
    motif_subfamilies = recognizers_by_level[4].subfamilies_by_uniprot_id(uniprot)
    comments = []

    [
      motif,
      pcm.length,
      pcm.consensus_string,
      uniprot,
      uniprot_infos.flat_map(&:uniprot_ac).join('; '),
      uniprot_infos.flat_map(&:primary_gene_name).join('; '),
      arity,
      quality,
      rank,
      num_words,
      best_auc_human, best_auc_mouse,
      num_datasets_human, num_datasets_mouse,
      release, source,
      motif_families.join(':separator:'),
      motif_subfamilies.join(':separator:'),
      uniprot_infos.flat_map(&:hgnc_ids).join('; '),
      uniprot_infos.flat_map(&:mgi_ids).join('; '),
      uniprot_infos.flat_map(&:entrezgene_ids).join('; '),
      comments.join(" "),
    ]
  }

  headers = [
    'Model name', 'Model length', 'Consensus', 'UniprotID', 'UniprotAC', 'Gene name',
    'Model type', 'Model quality', 'Model rank', 'Number of words in alignment',
    'Best AUC for human datasets', 'Best AUC for mouse datasets', 'Number of human datasets', 'Number of mouse datasets',
    'Release version', 'Data source type', 'Motif family', 'Motif subfamily',
    'HGNC', 'MGI', 'EntrezGene',
    'Comment'
  ]

  output_stream.puts headers.join("\t")
  model_infos.each do |infos|
    output_stream.puts infos.join("\t")
  end
end

def make_collection_summary(folder, species, arity, output_file)
  infos_by_uniprot_id = UniprotInfo
                        .each_in_file('uniprot_HomoSapiens_and_MusMusculus_lots_of_infos.tsv')
                        .group_by(&:uniprot_id)

 
  # we should get data for both species in order to know origin of motif obtained as cross-species
  original_motifs_origin = ['HUMAN', 'MOUSE'].flat_map{|species|
    File.readlines("hocomoco10/#{species}/#{arity}/final_collection.tsv").drop(1).map{|line|
      motif, hocomoco10_motif_origin = line.chomp.split("\t").values_at(0, 12)
      [motif, origin_by_motif_in_hocomoco10(hocomoco10_motif_origin)]
    }
  }.to_h

  motif_origin_infos = File.readlines("final_collection/#{arity}.tsv").map{|l|
    l.chomp.split("\t")[0,3]
  }.select{|origin, motif, original_motif|
    motif.match(/.*_#{species}\./)
  }
  original_motifs = motif_origin_infos.map{|origin, motif, original_motif| [motif, original_motif] }.to_h

  motif_names = Dir.glob("#{folder}/pcm/*").map{|fn|
    File.basename(fn, File.extname(fn))
  }.sort

  num_words_by_motif = motif_names.map{|motif|
    num_words = File.readlines("#{folder}/words/#{motif}.words").map(&:strip).reject(&:empty?).size
    [motif, num_words]
  }.to_h

  model_kind = ModelKind.get(arity)

  pcm_by_motif = motif_names.map{|motif|
    pcm = model_kind.read_pcm("#{folder}/pcm/#{motif}.#{model_kind.pcm_extension}")
    [motif, pcm]
  }.to_h

  File.open(output_file, 'w') do |fw|
    store_collection_summary(
      species, arity, motif_names,
      num_words_by_motif: num_words_by_motif,
      pcm_by_motif: pcm_by_motif,
      infos_by_uniprot_id: infos_by_uniprot_id,
      original_motifs: original_motifs,
      original_motifs_origin: original_motifs_origin,
      output_stream: fw
    )
  end
end


def repack_collection(species, arity, folder, motifs, hocomoco_prefix)
  model_kind = ModelKind.get(arity)
  rm_rf folder
  motifs.each do |motif|
    FileUtils.cp("final_collection/#{arity}/pcm/#{motif}.#{model_kind.pcm_extension}", "#{folder}/pcm/")
    FileUtils.cp("final_collection/#{arity}/pwm/#{motif}.#{model_kind.pwm_extension}", "#{folder}/pwm/")
    FileUtils.cp("final_collection/#{arity}/words/#{motif}.words", "#{folder}/words/")
    FileUtils.cp("final_collection/#{arity}/thresholds/#{motif}.thr", "#{folder}/thresholds/")
    FileUtils.cp(['direct', 'revcomp'].map{|orient| "final_collection/#{arity}/logo/#{motif}_#{orient}.png" }, "#{folder}/logo/")
    FileUtils.cp(['direct', 'revcomp'].map{|orient| "final_collection/#{arity}/logo_large/#{motif}_#{orient}.png" }, "#{folder}/logo_large/")
    FileUtils.cp(['direct', 'revcomp'].map{|orient| "final_collection/#{arity}/logo_small/#{motif}_#{orient}.png" }, "#{folder}/logo_small/")
  end
  make_collection_summary(folder, species, arity, File.join(folder, "#{hocomoco_prefix}final_collection_#{species}_#{arity}.tsv"))
  requested_pvalues = [0.001, 0.0005, 0.0001]
  thresholds_by_model = load_thresholds_by_model(folder, species, arity, requested_pvalues)
  save_standard_thresholds!(File.join(folder, "#{hocomoco_prefix}standard_thresholds_#{species}_#{arity}.txt"), thresholds_by_model, requested_pvalues)

  save_collection_in_single_files!(folder, species, arity, requested_pvalues, thresholds_by_model, hocomoco_prefix)
end

def archive_results(species, arity, folder, hocomoco_prefix)
  {
    "#{hocomoco_prefix}pcm_#{species}_#{arity}.tar.gz") => 'pcm',
    "#{hocomoco_prefix}pwm_#{species}_#{arity}.tar.gz") => 'pwm',
    "#{hocomoco_prefix}words_#{species}_#{arity}.tar.gz") => 'words',
    "#{hocomoco_prefix}thresholds_#{species}_#{arity}.tar.gz") => 'thresholds',
    "#{hocomoco_prefix}logo_#{species}_#{arity}.tar.gz") => 'logo',
    "#{hocomoco_prefix}logo_large_#{species}_#{arity}.tar.gz") => 'logo_large',
    "#{hocomoco_prefix}logo_small_#{species}_#{arity}.tar.gz") => 'logo_small',
  }.each{|archive_name, folder_to_pack|
    next  if File.exist?(archive_name)
    next  if !File.exist?(folder_to_pack)
    sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, archive_name), folder_to_pack
  }
end

def motifs_by_bundle(arity, species)
  result = Hash.new{|h,k| h[k] = [] }
  Dir.glob('final_collection/#{arity}/json/*.json').lazy.map{|fn|
    JSON.parse(File.read(fn), symbolize_names: true)
  }.select{|infos|
    infos[:species] == species
  }.each{|infos|
    infos[:bundle_list].each{|bundle|
      result[bundle] << infos[:name]
    }
  }
  result
end

desc 'Make final bundles'
task :repack_final_collection do
  ['HUMAN', 'MOUSE'].each do |species|
    ['mono', 'di'].each do |arity|
      motifs_by_bundle(arity, species).each{|bundle, motifs|
        repack_collection(species, arity, "final_bundle/hocomoco11/#{bundle}/#{species}/#{arity}", motifs, "HOCOMOCOv11_#{bundle}_")
      }
    end
  end
end

desc 'Archive final bundle subfolders'
task :archive_final_bundle do
  ['HUMAN', 'MOUSE'].each do |species|
    ['mono', 'di'].each do |arity|
      motifs_by_bundle(arity, species).each_key{|bundle|
        archive_results(species, arity, "final_bundle/hocomoco11/#{bundle}/#{species}/#{arity}", "HOCOMOCOv11_#{bundle}_")
      }
    end
  end
end
