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
def save_collection_in_single_files!(folder, species, arity, model_infos, requested_pvalues, thresholds_by_model)
  File.write File.join(folder, "pcms_#{species}_#{arity}.txt"), model_infos.map(&:pcm).map(&:to_s).join("\n")
  File.write File.join(folder, "pwms_#{species}_#{arity}.txt"), model_infos.map(&:pwm).map(&:to_s).join("\n")

  if arity == 'mono'
    File.write File.join(folder, "HOCOMOCOv10_#{species}_mono_meme_format.meme"), in_meme_format(model_infos.map(&:pcm))
    File.write File.join(folder, "HOCOMOCOv10_#{species}_mono_transfac_format.txt"), in_transfac_format(model_infos.map(&:pcm))
    File.write File.join(folder, "HOCOMOCOv10_#{species}_mono_jaspar_format.txt"), in_jaspar_format(model_infos.map(&:pcm))
    requested_pvalues.each do |requested_pvalue|
      File.write File.join(folder, "HOCOMOCOv10_#{species}_mono_homer_format_#{requested_pvalue}.motif"), in_homer_format(model_infos.map(&:pcm), thresholds_by_model, pvalue: requested_pvalue)
    end
  end
end

def calculate_thresholds_by_model(folder, model_infos, requested_pvalues)
  model_infos.map{|model_info|
    threshold_by_pvalue = Ape.run_find_threshold(
      File.join(folder, 'pwm', "#{model_info.full_name}.#{model_info.pwm_extension}"),
      requested_pvalues,
      discretization: 1000,
      mode: model_info.arity_type
    )
    [model_info.full_name, threshold_by_pvalue]
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

def calculate_all_thresholds!(folder, arity)
  sh 'java', '-cp', 'ape.jar',
      (arity == 'mono') ? 'ru.autosome.ape.PrecalculateThresholds' : 'ru.autosome.ape.di.PrecalculateThresholds',
      File.join(folder, "pwm"), File.join(folder, "thresholds"),
      '--pvalues', *['1e-15', '1.0', '1.01', 'mul'].join(','),
      '--discretization', 1000.to_s,
      '--silent'
end

def save_collection!(folder, model_infos, species, arity, requested_pvalues)
  mkdir_p File.join(folder, 'pcm')
  mkdir_p File.join(folder, 'pwm')
  mkdir_p File.join(folder, 'logo')
  mkdir_p File.join(folder, 'words')
  mkdir_p File.join(folder, 'logo_small')
  mkdir_p File.join(folder, 'logo_large')

  File.open(File.join(folder, "final_collection.html"), 'w') do |fw|
    print_html_table_by_model_infos(model_infos, stream: fw)
  end

  File.open(File.join(folder, "final_collection.tsv"), 'w') do |fw|
    print_csv_table(model_infos, stream: fw)
  end

  model_infos.each do |model_info|
    model_info.save_model_pack_into_folder!(folder)
  end

  thresholds_by_model = calculate_thresholds_by_model(folder, model_infos, requested_pvalues)
  save_standard_thresholds!(File.join(folder, "standard_thresholds_#{species}_#{arity}.txt"), thresholds_by_model, requested_pvalues)
  calculate_all_thresholds!(folder, arity)

  save_collection_in_single_files!(folder, species, arity, model_infos, requested_pvalues, thresholds_by_model)
end


desc 'Collect final collection'
task :make_final_collection do
  auc_infos_for_uniprot = AUCs.load_auc_infos_for_uniprot(min_weight_for_dataset: 0.65,
                                                          min_auc_for_model: 0.65)

  secondary_models = File.readlines('secondary_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }
  banned_models = File.readlines('banned_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }
  retracted_factors = File.readlines('retracted_factors.txt').map(&:chomp)

  best_models = collect_best_models(auc_infos_for_uniprot,
                                    secondary_models: secondary_models,
                                    banned_models: banned_models)

  quality_assessor = QualityAssessor.new(auc_infos_for_uniprot, best_models: best_models, secondary_models: secondary_models)

  File.open('final_collection.html', 'w') do |fw|
    print_html_table_for_grouped_models(
      auc_infos_for_uniprot,
      best_models.group_by(&:uniprot),
      quality_assessor,
      stream: fw
    )
  end


  to_be_reversed = File.readlines('revcomp_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }.to_set

  requested_pvalues = [0.001, 0.0005, 0.0001]

  rm_rf 'final_bundle'
  ['HUMAN', 'MOUSE'].each do |species|
    ['mono', 'di'].each do |arity|
      models = best_models.select{|model|
        model.species == species
      }.select{|model|
        model.arity_type == arity
      }

      # combine same models for several TFs into joint models
      all_model_infos = JointModel.grouped_models_from_scratch(models, auc_infos_for_uniprot, quality_assessor, to_be_reversed)
      good_model_infos = all_model_infos.reject{|model_info| retracted_factors.include?(model_info.uniprot) }
      retracted_model_infos = all_model_infos.select{|model_info| retracted_factors.include?(model_info.uniprot) }

      folder = "final_bundle/#{species}/#{arity}"
      save_collection!(folder, good_model_infos, species, arity, requested_pvalues)
      save_collection!("final_bundle/retracted/#{species}/#{arity}", retracted_model_infos, species, arity, requested_pvalues)

      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "pcm_#{species}_#{arity}.tar.gz"), 'pcm'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "pwm_#{species}_#{arity}.tar.gz"), 'pwm'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "words_#{species}_#{arity}.tar.gz"), 'words'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "thresholds_#{species}_#{arity}.tar.gz"), 'thresholds'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "logo_#{species}_#{arity}.tar.gz"), 'logo'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "logo_large_#{species}_#{arity}.tar.gz"), 'logo_large'
      sh 'tar', '-zhc', '-C', folder, '-f', File.join(folder, "logo_small_#{species}_#{arity}.tar.gz"), 'logo_small'
    end
  end
  if !Dir.glob('final_bundle/retracted/**/*').select{|fn| File.file?(fn)}.empty?
    sh 'tar', '-zhc', '-C', "final_bundle", '-f', "final_bundle/retracted.tar.gz", 'retracted'
  end
end
