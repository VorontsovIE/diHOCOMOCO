require 'thresholds_bsearch'

BACKGROUND_BY_SPECIES = {
  'HUMAN' => '0.09774531292656502,0.05049224075299731,0.07019109895771408,0.07682178619511619,0.0727342790964817,0.05203614856201394,0.010180820713495882,0.07019109895771408,0.059669884332282236,0.042565262995142815,0.05203614856201394,0.05049224075299731,0.06469420084013656,0.059669884332282236,0.0727342790964817,0.09774531292656502',
  'MOUSE' => '0.09124954151587066,0.05327746891945427,0.07340655447309075,0.07380976720188166,0.07444027240460285,0.0522326724288473,0.008258817805366036,0.07340655447309075,0.06218694059369016,0.04063209300331165,0.0522326724288473,0.05327746891945427,0.06371242131832879,0.06218694059369016,0.07444027240460285,0.09124954151587066',
}

def calculate_thresholds_cmd(pwm_fn, output_folder, background, additional_options)
  [
   'java', '-cp', 'ape.jar', 'ru.autosome.ape.di.PrecalculateThresholds',
    pwm_fn, output_folder, '--single-motif',
    '--background', background,
    '--pvalues', *['1e-15', '1.0', '1.01', 'mul'].join(','),
    '--discretization', 10000.to_s,
    *additional_options,
    '--silent',
  ].shelljoin
end

desc 'Threshold precalcuation for final bundle'
task :precalc_thresholds_for_final_bundle do
  ['mono', 'di'].each do |arity|
    output_folder = "final_collection/#{arity}/thresholds"
    FileUtils.mkdir_p(output_folder)
    pwm_folder = "final_collection/#{arity}/pwm"
    model_kind = ModelKind.get(arity)
    additional_options = (arity == 'mono') ? ['--from-mono'] : []

    motifs = Dir.glob("#{pwm_folder}/*").map{|fn|
      File.basename(fn, File.extname(fn))
    }

    motifs.reject{|motif|
      File.exist?("#{output_folder}/#{motif}.thr")
    }.each{|motif|
      species = motif.split('.').first.split('_').last
      pwm_fn = "#{pwm_folder}/#{motif}.#{model_kind.pwm_extension}"
      cmd = calculate_thresholds_cmd(pwm_fn, output_folder, BACKGROUND_BY_SPECIES[species], additional_options)
      puts cmd
    }
  end
end

task :put_thresholds_to_json do
  requested_pvalues = [0.001, 0.0005, 0.0001]
  ['mono', 'di'].each do |arity|
    motifs = Dir.glob("final_collection/#{arity}/json_processed/*.json").sort.each{|json_fn|
      motif_infos = JSON.parse(File.read(json_fn), symbolize_names: true)
      motif = motif_infos[:name]
      motif_infos[:threshold_pvalue_list] = load_threshold_pvalue_list("final_collection/#{arity}/thresholds/#{motif}.thr")
      motif_infos[:standard_thresholds] = thresholds_by_pvalues_bsearch(motif_infos[:threshold_pvalue_list], requested_pvalues)
      File.write(json_fn, motif_infos.to_json)
    }
  end
end
