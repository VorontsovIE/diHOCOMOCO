class QualityAssessor
  attr_reader :auc_infos_for_uniprot, :secondary_models
  def initialize(auc_infos_for_uniprot, secondary_models)
    @auc_infos_for_uniprot = auc_infos_for_uniprot
    @secondary_models = secondary_models
  end

  def hocomoco_quality(model)
    Models.hocomoco_qualities[model.model_name]
  end

  def aucs_by_dataset(model)
    auc_infos = auc_infos_for_uniprot[model.uniprot]
    auc_infos && auc_infos.aucs_for_model(model)
  end

  def num_datasets_passing_auc(model, threshold_auc)
    aucs_by_dataset(model).count{|dataset, auc| auc >= threshold_auc }
  end

  # has no validation at all
  #  (there were no control datasets or all datasets failed quality check)
  # or model was rejected because failed 0.65 threshold for weighted AUC
  def not_validated?(model)
    auc_infos = auc_infos_for_uniprot[model.uniprot]
    !auc_infos || !auc_infos.weighted_auc(model)
  end

  def calculate_quality(model)
    return 'S'  if secondary_models.include?(model)

    collection_name = model.collection_short_name
    is_hocomoco_model = (collection_name == 'HL')
    hocomoco_quality = is_hocomoco_model ? hocomoco_quality(model) : nil
    is_chipseq_model = Models::ChipseqCollections.include?(collection_name)

    if not_validated?(model)
      if is_hocomoco_model
        return hocomoco_quality
      else
        return 'D'
      end
    end

    num_datasets_pass_optimal_auc = num_datasets_passing_auc(model, 0.7)
    num_datasets_pass_minimal_auc = num_datasets_passing_auc(model, 0.65)

    if num_datasets_pass_optimal_auc >= 2
      $stderr.puts "#{model} has quality `A` but in hocomoco it had quality `#{hocomoco_quality}`" if is_hocomoco_model && hocomoco_quality != 'A'
      'A'
    elsif num_datasets_pass_optimal_auc == 1
      if is_chipseq_model && num_datasets_pass_minimal_auc >= 2
        $stderr.puts "#{model} has quality `B` but in hocomoco it had quality `#{hocomoco_quality}`" if is_hocomoco_model && hocomoco_quality != 'B'
        'B'
      elsif !is_chipseq_model
        $stderr.puts "#{model} has quality `B` but in hocomoco it had quality `#{hocomoco_quality}`" if is_hocomoco_model && hocomoco_quality != 'B'
        'B'
      else # is_chipseq_model && num_datasets_pass_minimal_auc == 1
        $stderr.puts "#{model} has quality `C` but in hocomoco it had quality `#{hocomoco_quality}`" if is_hocomoco_model && hocomoco_quality != 'C'
        'C'
      end
    else # num_datasets_pass_optimal_auc == 0
      if num_datasets_pass_minimal_auc >= 2
        $stderr.puts "#{model} has quality `C` but in hocomoco it had quality `#{hocomoco_quality}`" if is_hocomoco_model && hocomoco_quality != 'C'
        'C'
      else # num_datasets_pass_minimal_auc < 2
        if is_hocomoco_model
          hocomoco_quality(model)
        else
          'D'
        end
      end
    end
  end
end
