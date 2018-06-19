require 'csv'

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
  desc 'Rename motif collection alignment words into standardized ones (works on a stage of final bundle assembly)'
  task :rename_words do
    ['mono', 'di'].each do |model_kind|
      FileUtils.mkdir_p "final_collection/#{model_kind}/words"
      File.readlines("final_collection/#{model_kind}.tsv").drop(1).map{|line|
        line.chomp.split("\t").first(3)
      }.each{|origin, motif_final_name, motif_original_name|
        original_name = motif_original_name.split('~').last  # name w/o 'uniprot~collection~' part

        output_fn = "final_collection/#{model_kind}/words/#{motif_final_name}.words"
        if original_name.match(/\.H10(MO|DI)\./)
          FileUtils.cp "models/words/#{model_kind}/hocomoco_legacy/#{original_name}.words", output_fn
        else
          ext = (model_kind == 'mono') ? 'list' : 'dlist'
          words = File.readlines("models/words/#{model_kind}/chipseq/#{original_name}.#{ext}").drop(1).map{|l|
            l.chomp.split("\t")[2]
          }
          File.write(output_fn, words.map{|x| "#{x}\n" }.join)
        end
      }
    end
  end
end
