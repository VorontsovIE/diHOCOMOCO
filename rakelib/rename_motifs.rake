require 'csv'

MotifToUniprotMapping = Struct.new(:collection, :motif_name, :uniprot_ids, :entrezgene_ids) do
  def self.from_row(row) # row is an Array
    collection, motif_name, uniprot_ids, entrezgene_ids = *row
    uniprot_ids ||= ''
    entrezgene_ids ||= ''
    motif_name.gsub!('/','') # one of motifs have this shit in the name end in motif-uniprot FANTOM table
    # reject empty(represented as "-") and _RABIT uniprots
    uniprot_ids_normalized = uniprot_ids.split(',').map(&:strip).select{|uniprot_id| uniprot_id.match /_(HUMAN|MOUSE)/ }
    entrezgene_ids_normalized = entrezgene_ids.split(',').map(&:strip).reject(&:empty?)
    self.new(collection, motif_name, uniprot_ids_normalized, entrezgene_ids_normalized)
  end

  def self.each_in_file(filename, &block)
    CSV.foreach(filename, col_sep: "\t").drop(1).map{|row|
      self.from_row(row)
    }.each(&block)
  end

  def self.load_motif_to_uniprot_by_collection(filename)
    self.each_in_file(filename).group_by(&:collection).map{|collection, mot2uniprots|
      motif_mapping = mot2uniprots
                        .map{|mot2uniprot|
                          [mot2uniprot.motif_name, mot2uniprot.uniprot_ids]
                        }.to_h
      motif_mapping.default_proc = ->(h,k){ h[k] = [] }
      [collection, motif_mapping]
    }.to_h
  end
end


# if a motif has several uniprots, it's cloned several times (one for each uniprot)
def rename_motifs(src_glob, dest_folder,
                  short_collection_id:,
                  conv_to_uniprot_ids: ->(motif_name){
                    [ motif_name[/^.+_(HUMAN|MOUSE)/] ]
                  })
  mkdir_p dest_folder  unless Dir.exist?(dest_folder)
  FileList[src_glob].each do |src|
    extname = File.extname(src)
    motif_name = File.basename(src, extname)

    uniprot_ids = conv_to_uniprot_ids.call(motif_name)
    uniprot_ids.map{|uniprot_id|
      cp src, File.join(dest_folder, "#{uniprot_id}~#{short_collection_id}~#{motif_name}#{extname}")
    }
  end
end

namespace :collect_and_normalize_data do
  desc 'Rename motif collections into standardized ones; For motifs with several uniprots make single-uniprot copies'
  task :rename_motifs do
    motif_to_uniprot_mapping = MotifToUniprotMapping.load_motif_to_uniprot_by_collection('FANTOM5_phase2_KNOWN_motifs_IDmapping.txt')

    rename_motifs 'models/pcm/mono/hocomoco_legacy/*.pcm', 'models/pcm/mono/all/', short_collection_id: 'HL', conv_to_uniprot_ids: ->(motif){ motif_to_uniprot_mapping['HOCOMOCO'][motif] }
    rename_motifs 'models/pcm/mono/homer/*.pcm', 'models/pcm/mono/all/', short_collection_id: 'HO', conv_to_uniprot_ids: ->(motif){ motif_to_uniprot_mapping['HOMER'][motif] }
    rename_motifs 'models/pcm/mono/swissregulon/*.pcm', 'models/pcm/mono/all/', short_collection_id: 'SR', conv_to_uniprot_ids: ->(motif){ motif_to_uniprot_mapping['SWISSREGULON'][motif] }
    rename_motifs 'models/pcm/mono/jaspar/*.pcm', 'models/pcm/mono/all/', short_collection_id: 'JA', conv_to_uniprot_ids: ->(motif){ motif_to_uniprot_mapping['JASPAR'][motif] }

    rename_motifs 'models/pcm/mono/selex_ftr/*.pcm', 'models/pcm/mono/all/', short_collection_id: 'SE'
    rename_motifs 'models/pcm/mono/selex_rebuilt_ftr/*.pcm', 'models/pcm/mono/all/', short_collection_id: 'SMF'
    rename_motifs 'models/pcm/mono/selex_rebuilt_sub/*.pcm', 'models/pcm/mono/all/', short_collection_id: 'SMS'
    rename_motifs 'models/pcm/mono/chipseq/*.pcm', 'models/pcm/mono/all/', short_collection_id: 'CM'

    rename_motifs 'models/pcm/di/selex_rebuilt_ftr/*.dpcm', 'models/pcm/di/all/', short_collection_id: 'SDF'
    rename_motifs 'models/pcm/di/selex_rebuilt_sub/*.dpcm', 'models/pcm/di/all/', short_collection_id: 'SDS'
    rename_motifs 'models/pcm/di/chipseq/*.dpcm', 'models/pcm/di/all/', short_collection_id: 'CD'

    #####################
    # rename_motifs 'models/words/mono/hocomoco_legacy/*.words', 'models/words/mono/all/', short_collection_id: 'HL', conv_to_uniprot_ids: ->(motif){ motif_to_uniprot_mapping['HOCOMOCO'][motif] }
    # rename_motifs 'models/words/mono/homer/*.words', 'models/words/mono/all/', short_collection_id: 'HO', conv_to_uniprot_ids: ->(motif){ motif_to_uniprot_mapping['HOMER'][motif] }
    # rename_motifs 'models/words/mono/swissregulon/*.words', 'models/words/mono/all/', short_collection_id: 'SR', conv_to_uniprot_ids: ->(motif){ motif_to_uniprot_mapping['SWISSREGULON'][motif] }
    # rename_motifs 'models/words/mono/jaspar/*.words', 'models/words/mono/all/', short_collection_id: 'JA', conv_to_uniprot_ids: ->(motif){ motif_to_uniprot_mapping['JASPAR'][motif] }

    # rename_motifs 'models/words/mono/selex_ftr/*.words', 'models/words/mono/all/', short_collection_id: 'SE'
    # rename_motifs 'models/words/mono/selex_rebuilt_ftr/*.words', 'models/words/mono/all/', short_collection_id: 'SMF'
    # rename_motifs 'models/words/mono/selex_rebuilt_sub/*.words', 'models/words/mono/all/', short_collection_id: 'SMS'
    # rename_motifs 'models/words/mono/chipseq/*.words', 'models/words/mono/all/', short_collection_id: 'CM'

    # rename_motifs 'models/words/di/selex_rebuilt_ftr/*.words', 'models/words/di/all/', short_collection_id: 'SDF'
    # rename_motifs 'models/words/di/selex_rebuilt_sub/*.words', 'models/words/di/all/', short_collection_id: 'SDS'
    # rename_motifs 'models/words/di/chipseq/*.words', 'models/words/di/all/', short_collection_id: 'CD'
  end
end
