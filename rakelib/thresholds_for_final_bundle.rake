BACKGROUND_BY_SPECIES = {
  'HUMAN' => '0.09774531292656502,0.05049224075299731,0.07019109895771408,0.07682178619511619,0.0727342790964817,0.05203614856201394,0.010180820713495882,0.07019109895771408,0.059669884332282236,0.042565262995142815,0.05203614856201394,0.05049224075299731,0.06469420084013656,0.059669884332282236,0.0727342790964817,0.09774531292656502',
  'MOUSE' => '0.09124954151587066,0.05327746891945427,0.07340655447309075,0.07380976720188166,0.07444027240460285,0.0522326724288473,0.008258817805366036,0.07340655447309075,0.06218694059369016,0.04063209300331165,0.0522326724288473,0.05327746891945427,0.06371242131832879,0.06218694059369016,0.07444027240460285,0.09124954151587066',
}

def calculate_all_thresholds(folder, species, arity)
  additional_options = (arity == 'mono') ? ['--from-mono'] : []
  Dir.glob(File.join(folder, 'pwm', '*')).each do |fn|
    motif_name = File.basename(pwm_fn, File.extname(pwm_fn))
    output_folder = "#{folder}/thresholds/"
    output_fn = File.join(output_folder, "#{motif_name}.thr")
    unless File.exist?(output_fn)
      puts calculate_thresholds_cmd(pwm_fn, output_folder, BACKGROUND_BY_SPECIES[species], additional_options)
    end
  end
end

def calculate_thresholds_cmd(pwm_fn, output_folder, background, additional_options)
  [
   'java', '-cp', 'ape.jar', 'ru.autosome.ape.di.PrecalculateThresholds',
    fn, output_folder, '--single-motif',
    '--background', background,
    '--pvalues', *['1e-15', '1.0', '1.01', 'mul'].join(','),
    '--discretization', 10000.to_s,
    *additional_options,
    '--silent',
  ].shelljoin
end

desc 'Threshold precalcuation for final bundle'
task :precalc_thresholds_for_final_bundle do
  ['HUMAN', 'MOUSE'].each do |species|
    ['mono', 'di'].each do |arity|
      calculate_all_thresholds("final_bundle/hocomoco11/full/#{species}/#{arity}", species, arity)
    end
  end
end

desc 'Copy thresholds calculated for full collection into core collection'
task :precalc_thresholds_for_final_bundle_core do
  ['HUMAN', 'MOUSE'].each do |species|
    ['mono', 'di'].each do |arity|
      core_motifs = Dir.glob("final_bundle/hocomoco11/core/#{species}/#{arity}/pwm/*").map{|fn|
        File.basename(fn, File.extname(fn))
      }
      core_motifs.each{|motif|
        from = "final_bundle/hocomoco11/full/#{species}/#{arity}/thresholds/#{motif}.thr"
        to = "final_bundle/hocomoco11/core/#{species}/#{arity}/thresholds/#{motif}.thr"
        FileUtils.cp from, to
      }
    end
  end
end
