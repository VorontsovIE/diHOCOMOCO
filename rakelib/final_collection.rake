require 'models'
require 'json'
require 'information_content'
require 'auc_infos'
require 'motif_family_recognizer'
require 'ape_find_threshold'
require 'formatters'
require 'thresholds_bsearch'

# whole collection in a single file (one for all PCMs, one for all PWMs etc)
def save_collection_in_single_files!(folder, species, arity, motif_infos, requested_pvalues, hocomoco_prefix)
  File.open("#{folder}/#{hocomoco_prefix}pcms_#{species}_#{arity}.txt", 'w') {|fw|
    infos_for_motifs.each{|motif_infos|
      fw.puts ">#{motif_infos[:name]}\n"
      motif_infos[:pcm].each{|pos|
        fw.puts pos.join("\t")
      }
    }
  }
  File.open("#{folder}/#{hocomoco_prefix}pwms_#{species}_#{arity}.txt", 'w') {|fw|
    infos_for_motifs.each{|motif_infos|
      fw.puts ">#{motif_infos[:name]}\n"
      motif_infos[:pwm].each{|pos|
        fw.puts pos.join("\t")
      }
    }
  }

  if arity == 'mono'
    File.write File.join(folder, "#{hocomoco_prefix}#{species}_mono_meme_format.meme"), in_meme_format(infos_for_motifs)
    File.write File.join(folder, "#{hocomoco_prefix}#{species}_mono_transfac_format.txt"), in_transfac_format(infos_for_motifs)
    File.write File.join(folder, "#{hocomoco_prefix}#{species}_mono_jaspar_format.txt"), in_jaspar_format(infos_for_motifs)
    requested_pvalues.each do |requested_pvalue|
      File.write File.join(folder, "#{hocomoco_prefix}#{species}_mono_homer_format_#{requested_pvalue}.motif"), in_homer_format(infos_for_motifs, pvalue: requested_pvalue)
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

def motif_summary(motif, arity)
  motif_infos = JSON.parse(File.read("final_collection/#{arity}/json/#{motif}.json"), symbolize_names: true)
  uniprot_infos = motif_infos[:uniprot_infos]
  [
    motif,
    motif_infos[:length],
    motif_infos[:consensus_string],
    motif_infos[:uniprot],
    uniprot_infos[:uniprot_acs].join('; '),
    uniprot_infos[:primary_gene_names].join('; '),
    motif_infos[:model_kind], # aka `arity`
    motif_infos[:quality],
    motif_infos[:motif_index], # aka `rank`
    motif_infos[:num_words],
    motif_infos[:best_auc_human], motif_infos[:best_auc_mouse],
    motif_infos[:num_datasets_human], motif_infos[:num_datasets_mouse],
    motif_infos[:release], motif_infos[:source],
    motif_infos[:tfclass][:motif_families].join(':separator:'),
    motif_infos[:tfclass][:motif_subfamilies].join(':separator:'),
    uniprot_infos[:hgnc_ids].join('; '),
    uniprot_infos[:mgi_ids].join('; '),
    uniprot_infos[:entrezgene_ids].join('; '),
    motif_infos[:comments],
  ]
end

def make_collection_summary(folder, species, arity, output_file)
  motif_names = Dir.glob("#{folder}/pcm/*").map{|fn|
    File.basename(fn, File.extname(fn))
  }.sort

  headers = [
    'Model name', 'Model length', 'Consensus', 'UniprotID', 'UniprotAC', 'Gene name',
    'Model type', 'Model quality', 'Model rank', 'Number of words in alignment',
    'Best AUC for human datasets', 'Best AUC for mouse datasets', 'Number of human datasets', 'Number of mouse datasets',
    'Release version', 'Data source type', 'Motif family', 'Motif subfamily',
    'HGNC', 'MGI', 'EntrezGene',
    'Comment'
  ]

  File.open(output_file, 'w') do |fw|
    fw.puts headers.join("\t")
    motif_names.each{|motif|
      infos = motif_summary(motif, arity)
      fw.puts infos.join("\t")
    }
  end
end

def repack_collection(species, arity, infos_for_motifs, folder, hocomoco_prefix)
  model_kind = ModelKind.get(arity)
  rm_rf folder
  infos_for_motifs.each do |motif_infos|
    motif = motif_infos[:name]
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
  save_standard_thresholds!(File.join(folder, "#{hocomoco_prefix}standard_thresholds_#{species}_#{arity}.txt"), infos_for_motifs, requested_pvalues)

  save_collection_in_single_files!(folder, species, arity, motif_infos, requested_pvalues, hocomoco_prefix)
end

def archive_results(species, arity, folder, hocomoco_prefix)
  {
    "#{hocomoco_prefix}pcm_#{species}_#{arity}.tar.gz" => 'pcm',
    "#{hocomoco_prefix}pwm_#{species}_#{arity}.tar.gz" => 'pwm',
    "#{hocomoco_prefix}words_#{species}_#{arity}.tar.gz" => 'words',
    "#{hocomoco_prefix}thresholds_#{species}_#{arity}.tar.gz" => 'thresholds',
    "#{hocomoco_prefix}logo_#{species}_#{arity}.tar.gz" => 'logo',
    "#{hocomoco_prefix}logo_large_#{species}_#{arity}.tar.gz" => 'logo_large',
    "#{hocomoco_prefix}logo_small_#{species}_#{arity}.tar.gz" => 'logo_small',
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
      result[bundle] << infos
    }
  }

  result.each_key{|bundle|
    result[bundle].sort_by!{|motif_infos|
      motif_infos[:name]
    }
  }
  result
end

desc 'Make final bundles'
task :repack_final_collection do
  ['HUMAN', 'MOUSE'].each do |species|
    ['mono', 'di'].each do |arity|
      motifs_by_bundle(arity, species).each{|bundle, infos_for_motifs|
        repack_collection(species, arity, infos_for_motifs, "final_bundle/hocomoco11/#{bundle}/#{species}/#{arity}", "HOCOMOCOv11_#{bundle}_")
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
