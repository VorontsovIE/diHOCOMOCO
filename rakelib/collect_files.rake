def copy_files(src_glob, mapping)
  FileList[src_glob].each do |src|
    dest = src.pathmap(mapping)
    dir = File.dirname(dest)
    mkdir_p dir  unless Dir.exist?(dir)
    cp src, dest
  end
end

namespace :collect_and_normalize_data do
  desc 'Put models from different collections into standardified paths'
  task :collect_pcm do
    copy_files 'standard_motif_collections/pcm/hocomoco/*.pcm', 'models/pcm/mono/hocomoco_legacy/%n.pcm'
    copy_files 'standard_motif_collections/pcm/homer/*.pcm', 'models/pcm/mono/homer/%n.pcm'
    copy_files 'standard_motif_collections/pcm/swissregulon/*.pcm', 'models/pcm/mono/swissregulon/%n.pcm'
    copy_files 'standard_motif_collections/pcm/jaspar/*.pcm', 'models/pcm/mono/jaspar/%n.pcm'

    copy_files 'htselex/TAIPALE_uniprot.pcm/*.pat', 'models/pcm/mono/selex_ftr/%n.pcm'
    copy_files 'htselex_mono_di/mono_ftr/*.pcm', 'models/pcm/mono/selex_rebuilt_ftr/%n.pcm'
    copy_files 'htselex_mono_di/mono_sub/*.pcm', 'models/pcm/mono/selex_rebuilt_sub/%n.pcm'
    copy_files 'chipseq_models/mono/*.pcm', 'models/pcm/mono/chipseq/%n.pcm'

    copy_files 'chipseq_models/di/*.dpcm', 'models/pcm/di/chipseq/%n.dpcm'
    copy_files 'htselex_mono_di/di_ftr/*.dpcm', 'models/pcm/di/selex_rebuilt_ftr/%n.dpcm'
    copy_files 'htselex_mono_di/di_sub/*.dpcm', 'models/pcm/di/selex_rebuilt_sub/%n.dpcm'
  end
end
