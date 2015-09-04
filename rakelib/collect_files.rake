require 'jaspar'
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
  task :collect_pcm => :collect_pcm_jaspar do
    copy_files 'standard_motif_collections/pcm/hocomoco/*.pcm', 'models/pcm/mono/hocomoco_legacy/%n.pcm'
    copy_files 'standard_motif_collections/pcm/homer/*.pcm', 'models/pcm/mono/homer/%n.pcm'
    copy_files 'standard_motif_collections/pcm/swissregulon/*.pcm', 'models/pcm/mono/swissregulon/%n.pcm'
    # copy_files 'standard_motif_collections/pcm/jaspar/*.pcm', 'models/pcm/mono/jaspar/%n.pcm'

    copy_files 'pluripotency_factors/mono/*.pcm', 'models/pcm/mono/papatsenko/%n.pcm'
    copy_files 'pluripotency_factors/di/*.dpcm', 'models/pcm/di/papatsenko/%n.dpcm'

    copy_files 'htselex/TAIPALE_uniprot.pcm/*.pat', 'models/pcm/mono/selex_ftr/%n.pcm'
    copy_files 'htselex_mono_di/mono_ftr/*.pcm', 'models/pcm/mono/selex_rebuilt_ftr/%n.pcm'
    copy_files 'htselex_mono_di/mono_sub/*.pcm', 'models/pcm/mono/selex_rebuilt_sub/%n.pcm'
    copy_files 'selex_integrated/mono_fp/*.pcm', 'models/pcm/mono/selex_integrated/%n_HUMAN.pcm'
    copy_files 'chipseq_models/mono/*.pcm', 'models/pcm/mono/chipseq/%n.pcm'

    copy_files 'chipseq_models/di/*.dpcm', 'models/pcm/di/chipseq/%n.dpcm'
    copy_files 'htselex_mono_di/di_ftr/*.dpcm', 'models/pcm/di/selex_rebuilt_ftr/%n.dpcm'
    copy_files 'htselex_mono_di/di_sub/*.dpcm', 'models/pcm/di/selex_rebuilt_sub/%n.dpcm'
    copy_files 'selex_integrated/di_fp/*.dpcm', 'models/pcm/di/selex_integrated/%n_HUMAN.dpcm'
  end

  task :collect_pcm_jaspar do
    mkdir_p 'models/pcm/mono/jaspar/'
    uniprot_infos = UniprotInfo.each_in_file('uniprot_HomoSapiens_and_MusMusculus.txt').to_a
    jaspar_infos = Jaspar::Infos.new(
      position_counts_filename: 'standard_motif_collections_update/jaspar_2014/MATRIX_DATA.txt',
      taxonomy_filename: 'standard_motif_collections_update/jaspar_2014/TAX.txt',
      matrix_species_filename: 'standard_motif_collections_update/jaspar_2014/MATRIX_SPECIES.txt',
      matrix_proteins_filename: 'standard_motif_collections_update/jaspar_2014/MATRIX_PROTEIN.txt',
      matrix_name_infos_filename: 'standard_motif_collections_update/jaspar_2014/MATRIX.txt',
      uniprot_infos: uniprot_infos
    )
    jaspar_matrices = jaspar_infos.each_matrix.select{|info|
      # There is the only motif for Vertebrata taxon: HNF1A_HUMAN and is not linked to HNF1A_MOUSE
      # But there are no more HNF1A_MOUSE motifs so we don't care
      # ...and nevertheless it has no uniprot_id (because it can't be linked via provided identifiers)
      info.fit_species?('Homo sapiens') || info.fit_species?('Mus musculus') || info.fit_species?('Vertebrata')
    }.select{|info|
      info.collection == 'CORE'
    }.reject{|info|
      info.uniprot_ids.empty?
    }.group_by(&:name).each_value.map{|infos|
      infos.max_by(&:version)
    }.each{|info|
      info.uniprot_ids.each do |uniprot_id|
        mkdir_p File.join('models/pcm/mono/all/', uniprot_id)
        pcm_filename_1 = File.join('models/pcm/mono/jaspar/', "#{info.full_name}.pcm")
        pcm_filename_2 = File.join('models/pcm/mono/all/', uniprot_id, "#{uniprot_id}~JA~#{info.full_name}.pcm")
        $stderr.puts "Writing #{uniprot_id}~JA~#{info.full_name}"
        File.write(pcm_filename_1, info.matrix_str)
        File.write(pcm_filename_2, info.matrix_str(matrix_name: "#{uniprot_id}~JA~#{info.full_name}"))
      end
    }
  end
end
