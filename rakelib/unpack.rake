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
  sh 'tar', '-zxf', 'standard_motif_collections_update.tar.gz'
  ['pcm', 'pwm', 'ppm'].each do |model_type|
    mkdir_p File.join('standard_motif_collections', model_type)
    Dir.glob(File.join('standard_motif_collections', "*_#{model_type}.tar.gz")).each{|subarchive|
      sh 'tar', '-zxf', subarchive, '-C', File.join('standard_motif_collections', model_type)
    }
  end
end

desc 'Unpack archives'
task :unpack => [:unpack_zips, :unpack_motif_collections]
