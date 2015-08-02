module ModelKind
  def self.get(mode)
    case mode
    when Mono, Di
      mode
    when /^mono$/i
      @mono ||= Mono.new
    when /^di$/i
      @di ||= Di.new
    else
      raise "Unknown mode `#{mode.inspect}`; should be `:mono` or `:di`"
    end
  end

  class Mono
    def pwm_extension; 'pwm'; end
    def pcm_extension; 'pcm'; end
    def to_s; 'mono'; end
  end

  class Di
    def pwm_extension; 'dpwm'; end
    def pcm_extension; 'dpcm'; end
    def to_s; 'di'; end
  end
end

class Model
  attr_reader :uniprot, :collection_short_name, :model_name

  # AEBP2_HUMAN~CD~AEBP2_HUMAN^PEAKS030225, :di
  def initialize(model_fullname, mono_or_di_mode)
    @uniprot, @collection_short_name, @model_name = model_fullname.split('~')
    @mono_or_di_mode = ModelKind.get(mono_or_di_mode)
  end

  def pwm_extension; @mono_or_di_mode.pwm_extension; end
  def pwm_extension; @mono_or_di_mode.pcm_extension;; end

  def full_name
    [@uniprot, @collection_short_name, @model_name].join('~')
  end

  def path_to_pwm
    File.join('models/pwm', @mono_or_di_mode.to_s, 'all', uniprot, "#{full_name}.pwm")
  end

  def to_s; "<#{@mono_or_di_mode}: #{full_name}>"; end
  def inspect; to_s; end
end


module Models
  def self.mono_models
    @mono_models ||= FileList["models/pcm/mono/all/*/*.pcm"].pathmap('%n').map{|name| Model.new(name, :mono) }
  end

  def self.di_models
    @di_models ||= FileList["models/pcm/di/all/*/*.dpcm"].pathmap('%n').map{|name| Model.new(name, :di) }
  end

  def self.mono_models_by_uniprot(uniprot)
    @mono_models_by_uniprot ||= mono_models.group_by(&:uniprot)
    @mono_models_by_uniprot[uniprot]
  end

  def self.di_models_by_uniprot(uniprot)
    @di_models_by_uniprot ||= di_models.group_by(&:uniprot)
    @di_models_by_uniprot[uniprot]
  end

  def self.mono_uniprots
    mono_models.map(&:uniprot).uniq
  end

  def self.di_uniprots
    di_models.map(&:uniprot).uniq
  end
end
