def same_by?(models, &block)
  characteristics = models.map(&block)
  characteristics.all?{|ch| ch == characteristics.first }
end

# Calculate a characteristic for model. Make sure that all models have the same value or raise.
def take_and_check_consistency(models, &block)
  raise 'Can\'t take characteristic for empty model list'  if models.empty?
  if same_by?(models, &block)
    block.call(models.first)
  else
    raise 'Inconsistent characteristics for joint models #{models.inspect}'
  end
end

# Models should be equal. Normally models have one model, but sometimes
#   equal models can be related to different TFs
JointModel = Struct.new(:origin_models, :representative_model, :quality, :auc, :comments, :good_strand) do
  def initialize(*args, &block)
    super(*args, &block)
    raise "Not consistent models: #{origin_models.join(', ')}"  unless consistent?
  end

  def self.create_from_scratch(origin_models, auc_infos_for_uniprot, quality_assessor, to_be_reversed)
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

    good_strand = !to_be_reversed.include?(representative_model)

    self.new(origin_models, representative_model, quality, auc, comments, good_strand)
  end

  def consistent?
    origin_models.include?(representative_model) && \
    same_by?(origin_models, &:model_name) && \
    same_by?(origin_models, &:collection_short_name) && \
    same_by?(origin_models, &:arity_type) && \
    same_by?(origin_models, &:species) && \
    same_by?(origin_models){|m| m.pcm.matrix }
  end

  def species; representative_model.species; end
  def uniprot; representative_model.uniprot; end
  def arity_type; representative_model.arity_type; end
  def model_kind; ModelKind.get(arity_type); end

  def pcm_extension; model_kind.pcm_extension; end
  def pwm_extension; model_kind.pwm_extension; end

  def bundle_name
    {'mono' => 'H10MO', 'di' => 'H10DI'}[arity_type]
  end

  def full_name
    "#{uniprot}~#{bundle_name}~#{quality}"
  end

  # All models from the same collection with the same original name refer to the same model
  # We will join these models into one
  def self.grouped_models_from_scratch(models, auc_infos_for_uniprot, quality_assessor, to_be_reversed)
    models.group_by{|model|
      [model.collection_short_name, model.model_name].join('~')
    }.map{|_original_model_name, model_group|
      JointModel.create_from_scratch(model_group, auc_infos_for_uniprot, quality_assessor, to_be_reversed)
    }.sort_by(&:full_name)
  end

  def pcm
    (good_strand ? representative_model.pcm : representative_model.pcm.revcomp).named(full_name)
  end

  def pwm
    (good_strand ? representative_model.pwm : representative_model.pwm.revcomp).named(full_name)
  end

  def save_model_pack_into_folder!(folder)
    if good_strand
      FileUtils.cp representative_model.path_to_logo_direct, File.join(folder, 'logo', "#{full_name}_direct.png")
      FileUtils.cp representative_model.path_to_logo_revcomp, File.join(folder, 'logo', "#{full_name}_revcomp.png")
    else
      FileUtils.cp representative_model.path_to_logo_direct, File.join(folder, 'logo', "#{full_name}_revcomp.png")
      FileUtils.cp representative_model.path_to_logo_revcomp, File.join(folder, 'logo', "#{full_name}_direct.png")
    end

    File.write File.join(folder, 'pcm', "#{full_name}.#{pcm_extension}"), pcm.to_s
    File.write File.join(folder, 'pwm', "#{full_name}.#{pwm_extension}"), pwm.to_s
  end
end
