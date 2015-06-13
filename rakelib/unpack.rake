MotifToUniprotMapping = Struct.new(:collection, :motif_name, :uniprot_ids, :entrezgene_ids) do
  def self.from_row(row) # row is an Array
    collection, motif_name, uniprot_ids, entrezgene_ids = *row
    uniprot_ids ||= ''
    entrezgene_ids ||= ''
    # reject empty(represented as "-") and _RABIT uniprots
    motif_name.gsub!('/','') # one of motifs have this shit in motif-uniprot FANTOM table
    uniprot_ids_normalized = uniprot_ids.split(',').map(&:strip).select{|uniprot_id| uniprot_id.match /_(HUMAN|MOUSE)/ }
    entrezgene_ids_normalized = entrezgene_ids.split(',').map(&:strip).reject(&:empty?)
    new(collection, motif_name, uniprot_ids_normalized, entrezgene_ids_normalized)
  end
end

def copy_files(src_glob, mapping)
  FileList[src_glob].each do |src|
    dest = src.pathmap(mapping)
    dir = File.dirname(dest)
    mkdir_p dir  unless Dir.exist?(dir)
    cp src, dest
  end
end

def rename_motifs(src_glob, dest_folder, 
                  short_collection_id:,
                  conv_to_uniprot_ids: ->(motif){
                    [ motif[/^.+_(HUMAN|MOUSE)/] ]
                  })
  mkdir_p dest_folder  unless Dir.exist?(dest_folder)
  FileList[src_glob].each do |src|
    extname = File.extname(src)
    raise unless ['.pcm', '.dpcm'].include?(extname)
    motif = File.basename(src, extname)

    uniprot_ids = conv_to_uniprot_ids.call(motif)
    uniprot_ids.map{|uniprot_id|
      cp src, File.join(dest_folder, "#{uniprot_id}~#{short_collection_id}~#{motif}#{extname}")
    }
  end
end


# хорошо бы в macro-perfectos-ape JAVA встроить
def pcm_to_pwm(from_file, to_file)
  parser = Bioinform::MatrixParser.new(fix_nucleotides_number: 4)
  infos = parser.parse(File.read(from_file))
  name = infos[:name] || File.basename(from_file, '.pcm')
  # p from_file
  pcm = Bioinform::MotifModel::PCM.new(infos[:matrix]).named(name)

  pseudocount_calc = ->(pcm){
    max_count = pcm.each_position.map{|pos| pos.inject(0.0, &:+) }.max
    Math.log([max_count, 2].max)
  }
  converter = Bioinform::ConversionAlgorithms::PCM2PWMConverter.new(pseudocount: pseudocount_calc)
  pwm = converter.convert(pcm)

  File.write to_file, pwm.to_s
end

# хорошо бы в macro-perfectos-ape JAVA встроить
def dipcm_to_dipwm(from_file, to_file)
  parser = Bioinform::MatrixParser.new(fix_nucleotides_number: 16)
  infos = parser.parse(File.read(from_file))
  name = infos[:name] || File.basename(from_file, '.dpcm')
  pcm = Bioinform::MotifModel::DiPCM.new(infos[:matrix]).named(name)
  pwm = pcm.to_pwm
  File.write to_file, pwm.to_s
end

task :unpack_zips do
  Dir.glob('*.zip').each do |archive|
    output_dir = archive.pathmap('%X')
    mkdir_p output_dir
    sh 'unzip', archive, '-d', output_dir
  end
end

task :unpack_motif_collections do
  mkdir_p 'standard_motif_collections'
  sh 'tar', '-zxf', 'standard_motif_collections.tar.gz', '-C', 'standard_motif_collections'
  ['pcm', 'pwm', 'ppm'].each do |model_type|
    mkdir_p File.join('standard_motif_collections', model_type)
    Dir.glob(File.join('standard_motif_collections', "*_#{model_type}.tar.gz")).each{|subarchive|
      sh 'tar', '-zxf', subarchive, '-C', File.join('standard_motif_collections', model_type)
    }
  end
end

desc 'Unpack archives'
task :unpack => [:unpack_zips, :unpack_motif_collections]

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

desc 'Rename motif collections into standardized ones; For motifs with several uniprots make single-uniprot copies'
task :rename_motifs do
  require 'csv'

  # should I rename FASTA too?
  # ...

  motif_to_uniprot_mapping = CSV.foreach('FANTOM5_phase2_KNOWN_motifs_IDmapping.txt', col_sep: "\t").drop(1).map{|row|
    MotifToUniprotMapping.from_row(row)
  }.group_by(&:collection).map{|coll, mot2uniprots|
    motif_mapping = mot2uniprots \
                      #.reject{|mot2uniprot| mot2uniprot.uniprot_ids.empty? }
                      .map{|mot2uniprot|
                        [mot2uniprot.motif_name, mot2uniprot.uniprot_ids]
                      }.to_h
    [coll, (Hash.new{|h,k| h[k] = [] }).merge(motif_mapping)]
  }.to_h

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
end

desc 'Convert PCM --> PWM'
task :convert_pcm_to_pwm do
  require 'bioinform'
  require_relative '../lib/dipm'

  mkdir_p 'models/pwm/mono/all'
  FileList['models/pcm/mono/all/*.pcm'].each{|fn|
    pcm_to_pwm(fn, fn.pathmap('models/pwm/mono/all/%n.pwm'))
  }
  mkdir_p 'models/pwm/di/all'
  FileList['models/pcm/di/all/*.dpcm'].each{|fn|
    dipcm_to_dipwm(fn, fn.pathmap('models/pwm/di/all/%n.dpwm'))
  }
end
