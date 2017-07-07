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


# Designed to work with pcm-s/dipcm-s
# It not only renames file but also renames a motif (first line of matrix)
# If a motif has several uniprots, it's cloned several times (one for each uniprot)
def rename_motifs(src_glob, dest_folder,
                  short_collection_id:,
                  conv_to_uniprot_ids: ->(motif_name){
                    [ motif_name[/^.+_(HUMAN|MOUSE)/] ]
                  })
  mkdir_p dest_folder  unless Dir.exist?(dest_folder)
  FileList[src_glob].each do |src|
    extname = File.extname(src)
    motif_name = File.basename(src, extname)

    motif_text = File.readlines(src)
    motif_text = motif_text.first.match(/^>?\s*[a-zA-Z]/) ? motif_text.drop(1).join : motif_text.join
    uniprot_ids = conv_to_uniprot_ids.call(motif_name)
    uniprot_ids.each{|uniprot_id|
      subfolder = File.join(dest_folder, uniprot_id)
      mkdir_p(subfolder)  unless Dir.exist?(subfolder)
      motif_full_name = "#{uniprot_id}~#{short_collection_id}~#{motif_name}"
      dest = File.join(subfolder, "#{motif_full_name}#{extname}")
      next  if File.exist?(dest)
      $stderr.puts "Rename #{src} --> #{dest}"
      File.write(dest, "> #{motif_full_name}\n#{motif_text}")
    }
  end
end


def rename_words(src_glob, dest_folder,
                short_collection_id:,
                conv_to_uniprot_ids: ->(motif_name){
                  [ motif_name[/^.+_(HUMAN|MOUSE)/] ]
                })
  mkdir_p dest_folder  unless Dir.exist?(dest_folder)
  FileList[src_glob].each do |src|
    motif_name = File.basename(src, '.words')

    uniprot_ids = conv_to_uniprot_ids.call(motif_name)
    uniprot_ids.each{|uniprot_id|
      subfolder = File.join(dest_folder, uniprot_id)
      mkdir_p(subfolder)  unless Dir.exist?(subfolder)
      motif_full_name = "#{uniprot_id}~#{short_collection_id}~#{motif_name}"
      dest = File.join(subfolder, "#{motif_full_name}.words")
      next  if File.exist?(dest)
      $stderr.puts "Rename #{src} --> #{dest}"
      FileUtils.cp(src, dest)
    }
  end
end

namespace :collect_and_normalize_data do
  desc 'Rename motif collections into standardized ones; For motifs with several uniprots make single-uniprot copies'
  task :rename_motifs do
    rename_motifs 'models/pcm/mono/hocomoco_legacy/*.pcm', 'models/pcm/mono/all/', short_collection_id: 'HL'
    rename_motifs 'models/pcm/mono/chipseq/*.pcm', 'models/pcm/mono/all/', short_collection_id: 'CM'

    rename_motifs 'models/pcm/di/hocomoco_legacy/*.dpcm', 'models/pcm/di/all/', short_collection_id: 'DIHL'
    rename_motifs 'models/pcm/di/chipseq/*.dpcm', 'models/pcm/di/all/', short_collection_id: 'CD'
  end
end

namespace :collect_and_normalize_data do
  desc 'Rename motif collection alignment words into standardized ones; For motifs with several uniprots make single-uniprot copies'
  task :rename_words do
    motif_to_uniprot_mapping = MotifToUniprotMapping.load_motif_to_uniprot_by_collection('FANTOM5_phase2_KNOWN_motifs_IDmapping.txt')

    hocomoco_motif_to_uniprot = File.readlines('HOCOMOCOv9_motifs2uniprot.txt').drop(1).map(&:strip).reject(&:empty?).map{|line|
      motif, human_uniprots, mouse_uniprots = line.split("\t")
      [motif, (human_uniprots||'').split(',') + (mouse_uniprots||'').split(',')]
    }.to_h

    FileList['models/words/mono/hocomoco_legacy/*.words'].each do |fn|
      dest_folder = 'models/words/mono/all/'
      mkdir_p dest_folder  unless Dir.exist?(dest_folder)

      motif_name = File.basename(fn, '.words').sub(/_alignment$/, '')
      next  unless hocomoco_motif_to_uniprot.has_key?(motif_name)

      uniprot_ids = hocomoco_motif_to_uniprot[motif_name]
      uniprot_ids.each{|uniprot_id|
        subfolder = File.join(dest_folder, uniprot_id)
        mkdir_p(subfolder)  unless Dir.exist?(subfolder)
        motif_full_name = "#{uniprot_id}~HL~#{motif_name}"
        dest = File.join(subfolder, "#{motif_full_name}.words")
        next  if File.exist?(dest)
        $stderr.puts "Rename #{fn} --> #{dest}"
        # FileUtils.cp(fn, dest)
        words = File.readlines(fn).reject{|line|
          line.start_with?('#')
        }.map{|line|
          line.chomp.split("\t").first
        }
        File.write(dest, words.join("\n"))
      }
    end

    rename_words 'models/words/mono/chipseq/*.words', 'models/words/mono/all/', short_collection_id: 'CM'
    rename_words 'models/words/di/chipseq/*.words', 'models/words/di/all/', short_collection_id: 'CD'
  end
end
