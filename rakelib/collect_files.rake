require 'uniprot_info'

# mapping should map a file into a file name, not a folder name!
def copy_files(src_glob, mapping)
  FileList[src_glob].each do |src|
    dest = src.pathmap(mapping)
    next  if File.exist?(dest)
    dir = File.dirname(dest)
    mkdir_p dir  unless Dir.exist?(dir)
    cp src, dest
  end
end

namespace :collect_and_normalize_data do
  desc 'Put models from different collections into standardified paths'
  task :collect_pcm do
    copy_files 'standard_motif_collections/pcm/hocomoco10_human/*.pcm', 'models/pcm/mono/hocomoco_legacy/%n.pcm'
    copy_files 'standard_motif_collections/pcm/hocomoco10_mouse/*.pcm', 'models/pcm/mono/hocomoco_legacy/%n.pcm'
    copy_files 'chipseq_models_11/mono/*.pcm', 'models/pcm/mono/chipseq/%n.pcm'
  end


  desc 'Put words from different collections into standardified paths'
  task :collect_words do
    copy_files 'chipseq_models/mono/*.words', 'models/words/mono/chipseq/%n.words'
    copy_files 'chipseq_models/di/*.words', 'models/words/di/chipseq/%n.words'
  end
end
