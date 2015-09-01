class QualityAssessor
  def initialize(filtering)
    @hocomoco_qualities = self.class.load_hocomoco_qualities('hocomoco_qualities.tsv')
    @filtering = filtering
    @survived_models_by_uniprot = filtering.model_names.group_by{|model_name|
      Model.get_uniprot(model_name)
    }
  end

  def self.load_hocomoco_qualities(filename)
    File.readlines(filename).map{|line|  line.chomp.split("\t")  }.to_h
  end

  def num_datasets(model_name)
    @filtering.dataset_names_for_model(model_name).size
  end

  def num_survived_models(model_name)
    @survived_models_by_uniprot[Model.get_uniprot(model_name)].size
  end

  def hocomoco_quality(model_name)
    @hocomoco_qualities[Model.get_original_model_name(model_name)]
  end

  def num_datasets_passing_auc(model_name, auc_threshold)
    @filtering.aucs_for_model(model_name).count{|auc| auc >= auc_threshold }
  end

  def calculate_quality(model_name)
    is_hocomoco_model = model_name.match(/~HL~/)
    is_chipseq_model = model_name.match(/~CM~|~CD~/)
    is_selex_rebuilt_model = model_name.match(/~SMF~|~SMS~|~SDF~|~SDS~/)

    num_datasets_pass_highquality_auc = num_datasets_passing_auc(model_name, 0.9)
    num_datasets_pass_optimal_auc = num_datasets_passing_auc(model_name, 0.7)
    num_datasets_pass_minimal_auc = num_datasets_passing_auc(model_name, 0.6)
    
    # raise 'Impossible: best model doesn\'t pass filters'  if num_datasets(model_name) < 1
    # return hocomoco_quality(model_name)  if is_hocomoco_model

    if num_datasets_pass_optimal_auc >= 2
      return (num_datasets_pass_highquality_auc >= 1) ? 'A+' : 'A'
    elsif num_datasets_pass_optimal_auc == 1
      if !is_chipseq_model
        'B'
      elsif is_chipseq_model && num_datasets_pass_minimal_auc >= 2
        'B'
      else # is_chipseq_model && num_datasets_pass_minimal_auc == 1 
        'C'
      end
    else # num_datasets_pass_optimal_auc == 0
      if num_datasets_pass_minimal_auc >= 2 
        'C'
      elsif num_datasets_pass_minimal_auc == 1
        'D'
      else # num_datasets_pass_minimal_auc == 0
        if is_hocomoco_model
          hocomoco_quality(model_name)
        elsif is_selex_rebuilt_model
          'E'
        else
          'N/A'
        end
      end
    end
  end
end
