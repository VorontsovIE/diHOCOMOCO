require 'models'
require 'json'
require 'information_content'
require 'auc_infos'
require 'motif_family_recognizer'
require 'ape_find_threshold'
require 'formatters'
require 'thresholds_bsearch'

def save_annotation(infos_for_motifs, species, output_fn)
  header = [
    'Model',
    'Transcription factor',
    'Model length',
    'Quality',
    'Model rank',
    'Consensus',
    'Model release',
    'Data source',
    'Best auROC (human)',
    'Best auROC (mouse)',
    'Peak sets in benchmark (human)',
    'Peak sets in benchmark (mouse)',
    'Aligned words',
    'TF family',
    'TF subfamily',
    ((species == 'HUMAN') ? 'HGNC': 'MGI'),
    'EntrezGene',
    'UniProt ID',
    'UniProt AC',
  ]
  File.open(output_fn, 'w') do |fw|
    fw.puts header.join("\t")
    infos_for_motifs.each do |motif_infos|
      row = []
      row << motif_infos[:name]
      row << motif_infos[:uniprot_infos][:gene_names].join('; ')
      row += motif_infos.values_at(:length, :quality, :motif_index, :consensus_string, :release, :source)
      row += motif_infos.values_at(:best_auc_human, :best_auc_mouse, :num_datasets_human, :num_datasets_mouse,  :num_words)
      row += motif_infos[:tfclass].values_at(:motif_families, :motif_subfamilies).map{|xs| xs.join('; ') }
      if motif_infos[:species] == 'HUMAN'
        row << motif_infos[:uniprot_infos][:hgnc_ids].join('; ')
      else
        row << motif_infos[:uniprot_infos][:mgi_ids].join('; ')
      end
      row << motif_infos[:uniprot_infos][:entrezgene_ids].join('; ')
      row << motif_infos[:uniprot]
      row << motif_infos[:uniprot_infos][:uniprot_acs].join('; ')

      fw.puts row.join("\t")
    end
  end
end


# whole collection in a single file (one for all PCMs, one for all PWMs etc)
def save_collection_in_single_files!(folder, species, arity, infos_for_motifs, requested_pvalues, hocomoco_prefix)
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

  save_annotation(infos_for_motifs, species, File.join(folder, "#{hocomoco_prefix}annotation_#{species}_#{arity}.tsv"))
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
  motif_infos = JSON.parse(File.read("final_collection/#{arity}/json_processed/#{motif}.json"), symbolize_names: true)
  uniprot_infos = motif_infos[:uniprot_infos]
  [
    motif,
    motif_infos[:length],
    motif_infos[:consensus_string],
    motif_infos[:uniprot],
    uniprot_infos[:uniprot_acs].join('; '),
    uniprot_infos[:gene_names].join('; '),
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
  ['pcm', 'pwm', 'words', 'thresholds', 'logo', 'logo_small', 'logo_large'].each{|subfolder|
    FileUtils.mkdir_p "#{folder}/#{subfolder}/"
  }
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

  save_collection_in_single_files!(folder, species, arity, infos_for_motifs, requested_pvalues, hocomoco_prefix)
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
    next  if File.exist?(File.join(folder, archive_name))
    next  if !File.exist?(File.join(folder, folder_to_pack))
    sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, archive_name), folder_to_pack
  }
end

def motifs_by_bundle(arity, species)
  result = Hash.new{|h,k| h[k] = [] }
  Dir.glob("final_collection/#{arity}/json_processed/*.json").lazy.map{|fn|
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
  sh "tar -zhc -C final_bundle/hocomoco11/ -f final_bundle/hocomoco11/retracted.tar.gz retracted"
end

desc 'Fix pcm/pwm orientation; add consensus string to JSON configs'
task :improve_json do
  ['mono', 'di'].each do |arity|
    FileUtils.mkdir_p("final_collection/#{arity}/json_processed")
    model_kind = ModelKind.get(arity)
    Dir.glob("final_collection/#{arity}/json/*.json").each{|fn|
      motif_infos = JSON.parse(File.read(fn), symbolize_names: true)

      pcm = model_kind.read_pcm(motif_infos[:pcm_fn])
      pwm = model_kind.read_pwm(motif_infos[:pwm_fn])

      if motif_infos[:should_reverse]
        pcm = pcm.revcomp
        pwm = pwm.revcomp
      end

      words = motif_words(motif_infos[:words_fn])

      motif_infos[:consensus_string] = pcm.consensus_string
      motif_infos[:num_words] = words.size
      motif_infos[:pcm] = pcm.matrix
      motif_infos[:pwm] = pwm.matrix
      motif_infos[:words] = words

      motif_infos.delete(:pcm_fn)
      motif_infos.delete(:pwm_fn)
      motif_infos.delete(:words_fn)
      motif_infos.delete(:should_reverse)

      File.write("final_collection/#{arity}/json_processed/#{File.basename(fn)}", motif_infos.to_json)
    }
  end
end

desc 'Collect all hocomoco-11 json files into a single file'
task :collect_json_bundle_hocomoco11 do
  File.open('final_collection/hocomoco11_bundle.jsonl', 'w'){|fw|
    Dir.glob("final_collection/*/json/*.json").each{|fn| fw.puts File.read(fn) }
  }
end

desc 'Collect hocomoco-10 tsv tables into json a file'
task :collect_json_bundle_hocomoco10 do
  hocomoco10_infos = []
  ['mono', 'di'].each{|arity|
    ['HUMAN', 'MOUSE'].each{|species|
      File.readlines("hocomoco10/#{species}/#{arity}/final_collection.tsv").drop(1).map{|line|
        motif, hocomoco10_motif_origin = line.chomp.split("\t", 19)

        name, length, consensus, uniprot_id, uniprot_ac, gene_name, model_kind, quality, \
          num_words, weighted_auc, best_auc, datasets, hocomoco10_motif_origin, \
          motif_family, motif_subfamily, hgnc, mgi, entrezgene, comment = line.chomp.split("\t")
        raise 'Arity mismatch'  unless model_kind == arity
        raise 'Species mismatch'  unless motif.split('.').first.split('_').last == species
        info = {
          name: name, length: length, consensus: consensus,
          uniprot_id: uniprot_id, uniprot_ac: uniprot_ac,
          gene_names: gene_name.split('; ').reject(&:empty?),
          species: species, model_kind: model_kind, quality: quality,
          num_words: Integer(num_words),
          weighted_auc: (!weighted_auc || weighted_auc.empty?) ? nil : Float(weighted_auc),
          best_auc:     (!best_auc     || best_auc.empty?)     ? nil : Float(best_auc),
          datasets: datasets.split(', ').reject(&:empty?),
          original_motifs: hocomoco10_motif_origin.split(', ').reject(&:empty?),
          motif_families:     motif_family.split(':separator:').reject(&:empty?),
          motif_subfamilies:  motif_subfamily.split(':separator:').reject(&:empty?),
          hgnc_ids:              (hgnc || '').split('; ').reject(&:empty?),
          mgi_ids:                (mgi || '').split('; ').reject(&:empty?),
          entrezgene_ids:  (entrezgene || '').split('; ').reject(&:empty?),
          comment: comment || '',
        }
        hocomoco10_infos << info
      }
    }
  }

  File.open('final_collection/hocomoco10_bundle.jsonl', 'w'){|fw|
    hocomoco10_infos.each{|info|
      fw.puts(info.to_json)
    }
  }
end
