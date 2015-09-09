# Models should be equal. Normally models have one model, but sometimes
#   equal models can be related to different TFs
JointModel = Struct.new(:origin_models, :representative_model, :quality, :auc, :comments) do
  def initialize(*args, &block)
    super(*args, &block)
    raise "Not consistent models: #{origin_models.join(', ')}"  unless consistent?
  end

  def self.create_from_scratch(origin_models, auc_infos_for_uniprot, quality_assessor)
    origin_models = origin_models.sort_by(&:full_name)
    representative_model = origin_models.first

    quality = take_and_check_consistency(origin_models){|model|
      quality_assessor.calculate_quality(model)
    }

    aucs = origin_models.map{|model|
      auc_infos_for_uniprot[model.uniprot].weighted_auc(model)
    }.compact
    auc = aucs.empty? ? nil : median(aucs)

    comments = []
    if origin_models.size > 1
      comments << "Model is applicable to several TFs or complex subunits: #{origin_models.map(&:uniprot).join(', ')}."
    end

    comments << 'Secondary motif.' if quality == 'S'
    comments << 'Methylated DNA binding.'  if representative_model.model_name.match /!METH/

    self.new(origin_models, representative_model, quality, auc, comments)
  end

  def consistent?
    origin_models.include?(representative_model) && \
    same_by?(origin_models, &:model_name) && \
    same_by?(origin_models, &:collection_short_name) && \
    same_by?(origin_models, &:arity_type) && \
    same_by?(origin_models, &:species) && \
    same_by?(origin_models){|m| m.pcm.matrix }
  end

  def uniprot; representative_model.uniprot; end
  def arity_type; representative_model.arity_type; end
  def species; representative_model.arity_type; end

  def bundle_name
    {'mono' => 'H10MO', 'di' => 'H10DI'}[arity_type]
  end

  def full_name
    "#{uniprot}~#{bundle_name}~#{quality}"
  end

  # All models from the same collection with the same original name refer to the same model
  # We will join these models into one
  def self.grouped_models_from_scratch(models, auc_infos_for_uniprot, quality_assessor)
    models.group_by{|model|
      [model.collection_short_name, model.model_name].join('~')
    }.map{|_original_model_name, model_group|
      JointModel.create_from_scratch(model_group, auc_infos_for_uniprot, quality_assessor)
    }
  end
end
