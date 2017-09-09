require 'models'
require 'best_models'
require 'joint_model'
require 'auc_infos'
require 'quality_assessor'
require 'html_table_output'
require 'motif_family_recognizer'
require 'ape_find_threshold'
require 'formatters'

# whole collection in a single file (one for all PCMs, one for all PWMs etc)
def save_collection_in_single_files!(folder, species, arity, requested_pvalues, thresholds_by_model)
  model_kind = ModelKind.get(arity)
  pcms = Dir.glob("#{folder}/pcm/*").sort.map{|fn| model_kind.read_pcm(fn) }
  pwms = Dir.glob("#{folder}/pwm/*").sort.map{|fn| model_kind.read_pwm(fn) }

  File.write File.join(folder, "HOCOMOCOv11_pcms_#{species}_#{arity}.txt"), pcms.map(&:to_s).join("\n")
  File.write File.join(folder, "HOCOMOCOv11_pwms_#{species}_#{arity}.txt"), pwms.map(&:to_s).join("\n")

  if arity == 'mono'
    File.write File.join(folder, "HOCOMOCOv11_#{species}_mono_meme_format.meme"), in_meme_format(pcms)
    File.write File.join(folder, "HOCOMOCOv11_#{species}_mono_transfac_format.txt"), in_transfac_format(pcms)
    File.write File.join(folder, "HOCOMOCOv11_#{species}_mono_jaspar_format.txt"), in_jaspar_format(pcms)
    requested_pvalues.each do |requested_pvalue|
      File.write File.join(folder, "HOCOMOCOv11_#{species}_mono_homer_format_#{requested_pvalue}.motif"), in_homer_format(pcms, thresholds_by_model, pvalue: requested_pvalue)
    end
  end
end

def calculate_thresholds_by_model(folder, species, arity, requested_pvalues)
  Dir.glob(File.join(folder, '*')).map{|fn|
    threshold_by_pvalue = Ape.run_find_threshold(
      fn, requested_pvalues,
      discretization: 10000,
      background: BACKGROUND_BY_SPECIES[species],
      mode: 'di',
      additional_options: (arity == 'mono') ? ['--from-mono'] : []
    )
    [File.basename(fn, File.extname(fn)), threshold_by_pvalue]
  }.map{|model, threshold_by_pvalue|
    rounded_thresholds = threshold_by_pvalue.map{|pvalue, threshold|
      [pvalue, threshold.round(6)]
    }.to_h
    [model, rounded_thresholds]
  }.to_h
end

def save_standard_thresholds!(filename, thresholds_by_model, requested_pvalues)
  header = ['# P-values', *requested_pvalues]
  matrix = thresholds_by_model.map{|name, thresholds|
    [name, *thresholds.values_at(*requested_pvalues)]
  }
  File.write(filename, [header, *matrix].map{|row| row.join("\t") }.join("\n"))
end

BACKGROUND_BY_SPECIES = {
  'HUMAN' => '0.09774531292656502,0.05049224075299731,0.07019109895771408,0.07682178619511619,0.0727342790964817,0.05203614856201394,0.010180820713495882,0.07019109895771408,0.059669884332282236,0.042565262995142815,0.05203614856201394,0.05049224075299731,0.06469420084013656,0.059669884332282236,0.0727342790964817,0.09774531292656502',
  'MOUSE' => '0.09124954151587066,0.05327746891945427,0.07340655447309075,0.07380976720188166,0.07444027240460285,0.0522326724288473,0.008258817805366036,0.07340655447309075,0.06218694059369016,0.04063209300331165,0.0522326724288473,0.05327746891945427,0.06371242131832879,0.06218694059369016,0.07444027240460285,0.09124954151587066',
}

def calculate_all_thresholds!(folder, species, arity)
  additional_options = (arity == 'mono') ? ['--from-mono'] : []
  sh 'java', '-cp', 'ape.jar',
      'ru.autosome.ape.di.PrecalculateThresholds',
      File.join(folder, "pwm"), File.join(folder, "thresholds"),
      '--background', BACKGROUND_BY_SPECIES[species],
      '--pvalues', *['1e-15', '1.0', '1.01', 'mul'].join(','),
      '--discretization', 10000.to_s,
      *additional_options,
      '--silent'
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

def copy_by_glob(glob_from, folder_to)
  FileUtils.mkdir_p folder_to
  Dir.glob(glob_from).each{|fn|
    FileUtils.cp(fn, folder_to)
  }
end

def get_auc_stats(auc_fn)
  return({ best_auc: nil, num_datasets: nil })  unless File.exist?(auc_fn)
  aucs = File.readlines(auc_fn).map{|l|
    Float(l.chomp.split("\t")[1])
  }
  { best_auc: aucs.max, num_datasets: aucs.size }
end

def get_release_and_source(original_motif)
  if original_motif.match(/~(CM|CD)~/)
    release = 'HOCOMOCOv11'
    source = 'ChIP-Seq'
  else
    prev_motif = original_motif.split('~').last
    if original_motifs_origin[prev_motif] == 'HOCOMOCO v9'
      release = 'HOCOMOCOv9'
      source = 'Integrative'
    else
      release = 'HOCOMOCOv10'
      source = original_motifs_origin[prev_motif]
    end
  end
end

desc 'Collect final collection'
task :repack_final_collection do
  requested_pvalues = [0.001, 0.0005, 0.0001]
  infos_by_uniprot_id = UniprotInfo
                        .each_in_file('uniprot_HomoSapiens_and_MusMusculus_lots_of_infos.tsv')
                        .group_by(&:uniprot_id)

  original_motifs_origin = ['HUMAN', 'MOUSE'].flat_map{|species|
    ['mono', 'di'].flat_map{|arity|
      File.readlines("hocomoco10/#{species}/#{arity}/final_collection.tsv").drop(1).map{|line|
        motif, hocomoco10_motif_origin = line.chomp.split("\t").values_at(0, 12)
        [motif, origin_by_motif_in_hocomoco10(hocomoco10_motif_origin)]
      }
    }
  }.to_h

  rm_rf 'final_bundle'
  ['HUMAN', 'MOUSE'].each do |species|
    ['mono', 'di'].each do |arity|
      folder = "final_bundle/#{species}/#{arity}"

      copy_by_glob("final_collection/#{arity}/pcm/*_#{species}.*", "#{folder}/pcm")
      copy_by_glob("final_collection/#{arity}/pwm/*_#{species}.*", "#{folder}/pwm")
      copy_by_glob("final_collection/#{arity}/words/*_#{species}.*", "#{folder}/words")
      copy_by_glob("final_collection/#{arity}/logo/*_#{species}.*", "#{folder}/logo")
      copy_by_glob("final_collection/#{arity}/logo_large/*_#{species}.*", "#{folder}/logo_large")
      copy_by_glob("final_collection/#{arity}/logo_small/*_#{species}.*", "#{folder}/logo_small")

      thresholds_by_model = calculate_thresholds_by_model("#{folder}/pwm", species, arity, requested_pvalues)
      save_standard_thresholds!(File.join(folder, "standard_thresholds_#{species}_#{arity}.txt"), thresholds_by_model, requested_pvalues)

      calculate_all_thresholds!(folder, species, arity)
      save_collection_in_single_files!(folder, species, arity, requested_pvalues, thresholds_by_model)

      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "pcm_#{species}_#{arity}.tar.gz"), 'pcm'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "pwm_#{species}_#{arity}.tar.gz"), 'pwm'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "words_#{species}_#{arity}.tar.gz"), 'words'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "thresholds_#{species}_#{arity}.tar.gz"), 'thresholds'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "logo_#{species}_#{arity}.tar.gz"), 'logo'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "logo_large_#{species}_#{arity}.tar.gz"), 'logo_large'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "logo_small_#{species}_#{arity}.tar.gz"), 'logo_small'

      motif_origin_infos = File.readlines("final_collection/#{arity}.tsv").map{|l|
        l.chomp.split("\t")[0,3]
      }.select{|origin, motif, original_motif|
        motif.match(/.*_#{species}\./)
      }
      original_motifs = motif_origin_infos.map{|origin, motif, original_motif| [motif, original_motif] }.to_h
      origins =  motif_origin_infos.map{|origin, motif, original_motif| [motif, origin] }.to_h

      recognizers_by_level = PROTEIN_FAMILY_RECOGNIZERS[species]

      model_kind = ModelKind.get(arity)
      motif_names = Dir.glob("final_bundle/#{species}/#{arity}/pcm/*").map{|fn|
        File.basename(fn, File.extname(fn))
      }.sort

      model_infos = motif_names.map{|motif|
        pcm = model_kind.read_pcm("final_collection/#{arity}/pcm/#{motif}.#{model_kind.pcm_extension}")
        num_words = File.readlines("#{folder}/words/#{motif}.words").map(&:strip).reject(&:empty?).size

        best_auc_human, num_datasets_human = get_auc_stats("auc/#{arity}/HUMAN_datasets/#{original_motifs[motif]}.txt").values_at(:best_auc, :num_datasets)
        best_auc_mouse, num_datasets_mouse = get_auc_stats("auc/#{arity}/MOUSE_datasets/#{original_motifs[motif]}.txt").values_at(:best_auc, :num_datasets)

        original_motif = original_motifs[motif]
        release, source = get_release_and_source(original_motif)

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


      File.open(File.join(folder, "final_collection.tsv"), 'w') do |fw|
        fw.puts headers.join("\t")
        model_infos.each do |infos|
          fw.puts infos.join("\t")
        end
      end
    end
  end
end
