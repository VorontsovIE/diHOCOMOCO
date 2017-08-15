require 'models'
require 'best_models'
require 'joint_model'
require 'auc_infos'
require 'quality_assessor'
require 'html_table_output'
require 'motif_family_recognizer'
require 'ape_find_threshold'
require 'formatters'

# whole collection in a single file (one for all PCMs, one for all PWMs etc)
def save_collection_in_single_files!(folder, species, arity, requested_pvalues, thresholds_by_model)
  model_kind = ModelKind.get(arity)
  pcms = Dir.glob("#{folder}/pcm/*").sort.map{|fn| model_kind.read_pcm(fn) }
  pwms = Dir.glob("#{folder}/pwm/*").sort.map{|fn| model_kind.read_pwm(fn) }

  File.write File.join(folder, "HOCOMOCOv11_pcms_#{species}_#{arity}.txt"), pcms.map(&:to_s).join("\n")
  File.write File.join(folder, "HOCOMOCOv11_pwms_#{species}_#{arity}.txt"), pwms.map(&:to_s).join("\n")

  if arity == 'mono'
    File.write File.join(folder, "HOCOMOCOv11_#{species}_mono_meme_format.meme"), in_meme_format(pcms)
    File.write File.join(folder, "HOCOMOCOv11_#{species}_mono_transfac_format.txt"), in_transfac_format(pcms)
    File.write File.join(folder, "HOCOMOCOv11_#{species}_mono_jaspar_format.txt"), in_jaspar_format(pcms)
    requested_pvalues.each do |requested_pvalue|
      File.write File.join(folder, "HOCOMOCOv11_#{species}_mono_homer_format_#{requested_pvalue}.motif"), in_homer_format(pcms, thresholds_by_model, pvalue: requested_pvalue)
    end
  end
end

def calculate_thresholds_by_model(folder, species, arity, requested_pvalues)
  Dir.glob(File.join(folder, '*')).map{|fn|
    threshold_by_pvalue = Ape.run_find_threshold(
      fn, requested_pvalues,
      discretization: 10000,
      background: BACKGROUND_BY_SPECIES[species],
      mode: 'di',
      additional_options: (arity == 'mono') ? ['--from-mono'] : []
    )
    [File.basename(fn, File.extname(fn)), threshold_by_pvalue]
  }.map{|model, threshold_by_pvalue|
    rounded_thresholds = threshold_by_pvalue.map{|pvalue, threshold|
      [pvalue, threshold.round(6)]
    }.to_h
    [model, rounded_thresholds]
  }.to_h
end

def save_standard_thresholds!(filename, thresholds_by_model, requested_pvalues)
  header = ['# P-values', *requested_pvalues]
  matrix = thresholds_by_model.map{|name, thresholds|
    [name, *thresholds.values_at(*requested_pvalues)]
  }
  File.write(filename, [header, *matrix].map{|row| row.join("\t") }.join("\n"))
end

BACKGROUND_BY_SPECIES = {
  'HUMAN' => '0.09774531292656502,0.05049224075299731,0.07019109895771408,0.07682178619511619,0.0727342790964817,0.05203614856201394,0.010180820713495882,0.07019109895771408,0.059669884332282236,0.042565262995142815,0.05203614856201394,0.05049224075299731,0.06469420084013656,0.059669884332282236,0.0727342790964817,0.09774531292656502',
  'MOUSE' => '0.09124954151587066,0.05327746891945427,0.07340655447309075,0.07380976720188166,0.07444027240460285,0.0522326724288473,0.008258817805366036,0.07340655447309075,0.06218694059369016,0.04063209300331165,0.0522326724288473,0.05327746891945427,0.06371242131832879,0.06218694059369016,0.07444027240460285,0.09124954151587066',
}

def calculate_all_thresholds!(folder, arity)
  additional_options = (arity == 'mono') ? ['--from-mono'] : []
  sh 'java', '-cp', 'ape.jar',
      'ru.autosome.ape.di.PrecalculateThresholds',
      File.join(folder, "pwm"), File.join(folder, "thresholds"),
      '--background', BACKGROUND_BY_SPECIES[species],
      '--pvalues', *['1e-15', '1.0', '1.01', 'mul'].join(','),
      '--discretization', 10000.to_s,
      *additional_options,
      '--silent'
end

desc 'Collect final collection'
task :repack_final_collection do
  requested_pvalues = [0.001, 0.0005, 0.0001]
  rm_rf 'final_bundle'
  ['HUMAN', 'MOUSE'].each do |species|
    ['mono', 'di'].each do |arity|
      folder = "final_bundle/#{species}/#{arity}"
      FileUtils.mkdir_p "#{folder}/pcm"
      FileUtils.mkdir_p "#{folder}/pwm"
      FileUtils.mkdir_p "#{folder}/words"
      FileUtils.mkdir_p "#{folder}/logo"
      FileUtils.mkdir_p "#{folder}/logo_large"
      FileUtils.mkdir_p "#{folder}/logo_small"
      Dir.glob("final_collection/#{arity}/pcm/*_#{species}.*").each{|fn| FileUtils.cp(fn, "#{folder}/pcm") }
      Dir.glob("final_collection/#{arity}/pwm/*_#{species}.*").each{|fn| FileUtils.cp(fn, "#{folder}/pwm") }
      Dir.glob("final_collection/#{arity}/words/*_#{species}.*").each{|fn| FileUtils.cp(fn, "#{folder}/words") }
      Dir.glob("final_collection/#{arity}/logo/*_#{species}.*").each{|fn| FileUtils.cp(fn, "#{folder}/logo") }
      Dir.glob("final_collection/#{arity}/logo_large/*_#{species}.*").each{|fn| FileUtils.cp(fn, "#{folder}/logo_large") }
      Dir.glob("final_collection/#{arity}/logo_small/*_#{species}.*").each{|fn| FileUtils.cp(fn, "#{folder}/logo_small") }

      thresholds_by_model = calculate_thresholds_by_model("#{folder}/pwm", species, arity, requested_pvalues)
      save_standard_thresholds!(File.join(folder, "standard_thresholds_#{species}_#{arity}.txt"), thresholds_by_model, requested_pvalues)

      calculate_all_thresholds!(folder, arity)
      save_collection_in_single_files!(folder, species, arity, requested_pvalues, thresholds_by_model)

      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "pcm_#{species}_#{arity}.tar.gz"), 'pcm'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "pwm_#{species}_#{arity}.tar.gz"), 'pwm'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "words_#{species}_#{arity}.tar.gz"), 'words'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "thresholds_#{species}_#{arity}.tar.gz"), 'thresholds'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "logo_#{species}_#{arity}.tar.gz"), 'logo'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "logo_large_#{species}_#{arity}.tar.gz"), 'logo_large'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "logo_small_#{species}_#{arity}.tar.gz"), 'logo_small'
    end
  end
end
