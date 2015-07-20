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

desc 'Move all motifs/words into their folders, normalize their names, convert to PWMs and so on'
task 'collect_and_normalize_data' do
  Rake::Task['collect_and_normalize_data:collect_pcm'].invoke
  Rake::Task['collect_and_normalize_data:rename_motifs'].invoke
  Rake::Task['collect_and_normalize_data:convert_pcm_to_pwm'].invoke
end
