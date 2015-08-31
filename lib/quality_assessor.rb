class QualityAssessor
  def initialize(filtering)
    @hocomoco_qualities = self.class.load_hocomoco_qualities('hocomoco_genes_infos.csv')
    @filtering = filtering
    @survived_models_by_uniprot = filtering.model_names.group_by{|model_name|
      Model.get_uniprot(model_name)
    }
  end

  def self.load_hocomoco_qualities(filename)
    File.readlines(filename).drop(1).map{|line|
      line.chomp.split("\t")
    }.map{|row|
      [row[0],row[2]]
    }.to_h
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

  def calculate_quality(model_name)
    is_hocomoco_model = model_name.match(/~HL~/)
    is_chipseq_model = model_name.match(/~CM~|~CD~/)
    
    raise 'Impossible: best model doesn\'t pass filters'  if num_datasets(model_name) < 1
    return hocomoco_quality(model_name)  if is_hocomoco_model

    if num_datasets(model_name) >= 2
      'A'
    else      
      if is_chipseq_model
        (num_survived_models(model_name) > 1) ? 'C' : 'D'
      else
        'B'
      end
    end
  end
end
