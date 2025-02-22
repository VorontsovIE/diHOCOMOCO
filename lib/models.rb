require 'bioinform'
require_relative 'dipm'
require_relative 'information_content'

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

  class Base
    def initialize; raise NotImplementedError, "abstract class can't be instantiated"; end
    def to_s; arity_type; end
    def inspect; arity_type; end
    def read_pcm(path_to_pcm, only_matrix: false)
      parser = Bioinform::MatrixParser.new(fix_nucleotides_number: alphabet_size, nucleotides_in: :columns)
      infos = parser.parse(File.read(path_to_pcm))
      name = infos[:name] || File.basename(path_to_pcm, ".#{pcm_extension}")
      matrix = create_pcm(infos[:matrix])
      only_matrix ? infos[:matrix] : create_pcm(infos[:matrix]).named(name)
    end
    def read_pwm(path_to_pwm, only_matrix: false)
      parser = Bioinform::MatrixParser.new(fix_nucleotides_number: alphabet_size, nucleotides_in: :columns)
      infos = parser.parse(File.read(path_to_pwm))
      name = infos[:name] || File.basename(path_to_pwm, ".#{pwm_extension}")
      only_matrix ? infos[:matrix] : create_pwm(infos[:matrix]).named(name)
    end
  end

  class Mono < Base
    def initialize; end
    def alphabet_size; 4; end
    def arity_type; 'mono'; end
    def pwm_extension; 'pwm'; end
    def pcm_extension; 'pcm'; end
    def wordlist_extension; 'list'; end
    def sarus_class; 'ru.autosome.SARUS'; end
    def create_pcm(matrix); Bioinform::MotifModel::PCM.new(matrix); end
    def create_pwm(matrix); Bioinform::MotifModel::PWM.new(matrix); end
  end

  class Di < Base
    def initialize; end
    def alphabet_size; 16; end
    def arity_type; 'di'; end
    def pwm_extension; 'dpwm'; end
    def pcm_extension; 'dpcm'; end
    def wordlist_extension; 'dlist'; end
    def sarus_class; 'ru.autosome.di.SARUS'; end
    def create_pcm(matrix); Bioinform::MotifModel::DiPCM.new(matrix); end
    def create_pwm(matrix); Bioinform::MotifModel::DiPWM.new(matrix); end
  end
end

class Model
  attr_reader :uniprot, :collection_short_name, :model_name, :mono_or_di_mode

  # AEBP2_HUMAN~CD~AEBP2_HUMAN^PEAKS030225, :di
  def initialize(model_fullname, mono_or_di_mode)
    @uniprot, @collection_short_name, @model_name = model_fullname.split('~')
    @mono_or_di_mode = ModelKind.get(mono_or_di_mode)
  end

  def ==(other)
    other.is_a?(Model) && \
    @uniprot == other.uniprot && \
    @collection_short_name == other.collection_short_name && \
    @model_name == other.model_name && \
    @mono_or_di_mode == other.mono_or_di_mode
  end
  def eql?(other)
    other.class == self.class && self == other
  end
  def hash
    [@uniprot, @collection_short_name, @model_name, @mono_or_di_mode].hash
  end

  def <=>(other)
    full_name <=> other.full_name
  end

  def pcm; @pcm ||= @mono_or_di_mode.read_pcm(path_to_pcm); end
  def pwm; @pwm ||= @mono_or_di_mode.read_pwm(path_to_pwm); end
  def length; pcm.length; end
  def arity_type; @mono_or_di_mode.arity_type; end

  def pwm_extension; @mono_or_di_mode.pwm_extension; end
  def pcm_extension; @mono_or_di_mode.pcm_extension; end

  def full_name
    [@uniprot, @collection_short_name, @model_name].join('~')
  end

  def semiuniprot; @uniprot.split('_').first; end # species-independent protein name
  def species; @uniprot.split('_').last; end

  def path_to_pcm
    File.join('models/pcm', @mono_or_di_mode.to_s, 'all', uniprot, "#{full_name}.#{pcm_extension}")
  end

  def path_to_pwm
    File.join('models/pwm', @mono_or_di_mode.to_s, 'all', uniprot, "#{full_name}.#{pwm_extension}")
  end

  def path_to_words
    File.join('models/words', @mono_or_di_mode.to_s, 'all', uniprot, "#{full_name}.words")
  end

  def num_words_in_alignment
    File.readlines(path_to_words).map(&:chomp).reject(&:empty?).size
  end

  def path_to_logo; path_to_logo_direct; end

  def path_to_logo_direct
    File.join('models/logo', uniprot, "#{full_name}_direct.png")
  end

  def path_to_logo_revcomp
    File.join('models/logo', uniprot, "#{full_name}_revcomp.png")
  end

  def path_to_small_logo_direct
    File.join('models/logo_small', uniprot, "#{full_name}_direct.png")
  end

  def path_to_small_logo_revcomp
    File.join('models/logo_small', uniprot, "#{full_name}_revcomp.png")
  end

  def path_to_large_logo_direct
    File.join('models/logo_large', uniprot, "#{full_name}_direct.png")
  end

  def path_to_large_logo_revcomp
    File.join('models/logo_large', uniprot, "#{full_name}_revcomp.png")
  end

  def to_s; "<#{@mono_or_di_mode}: #{full_name}>"; end
  def inspect; to_s; end

  def self.get_original_model_name(model_fullname)
    model_fullname.split('~')[2]
  end
  def self.get_uniprot(model_fullname)
    model_fullname.split('~')[0]
  end
  def self.get_collection_short_name(model_fullname)
    model_fullname.split('~')[1]
  end

  def self.new_by_name(fullname)
    collection = fullname.split('~')[1]
    in_mono_collection = Models::MonoCollections.include?(collection)
    in_di_collection = Models::DiCollections.include?(collection)
    raise "Oops. Collection `#{collection}` is both mono and dinucleotide" if in_mono_collection && in_di_collection
    raise "Oops. Collection `#{collection}` is neither mono nor dinucleotide" if !in_mono_collection && !in_di_collection
    self.new(fullname, (in_mono_collection ? :mono : :di))
  end
end


module Models
  def self.models_by_type(model_kind)
    glob = "models/pcm/#{model_kind.arity_type}/all/*/*.#{model_kind.pcm_extension}"
    FileList[glob].pathmap('%n').uniq.sort.map{|name|
      Model.new(name, model_kind.arity_type)
    }
  end

  def self.mono_models
    @mono_models ||= FileList["models/pcm/mono/all/*/*.pcm"].pathmap('%n').map{|name| Model.new(name, :mono) }
  end

  def self.di_models
    @di_models ||= FileList["models/pcm/di/all/*/*.dpcm"].pathmap('%n').map{|name| Model.new(name, :di) }
  end

  def self.all_models
    @all_models ||= mono_models + di_models
  end

  def self.mono_models_by_uniprot(uniprot)
    @mono_models_by_uniprot ||= mono_models.group_by(&:uniprot)
    @mono_models_by_uniprot.fetch(uniprot) { [] }
  end

  def self.di_models_by_uniprot(uniprot)
    @di_models_by_uniprot ||= di_models.group_by(&:uniprot)
    @di_models_by_uniprot.fetch(uniprot) { [] }
  end

  def self.all_models_by_uniprot(uniprot)
    @all_models_by_uniprot ||= all_models.group_by(&:uniprot)
    @all_models_by_uniprot.fetch(uniprot) { [] }
  end

  def self.mono_uniprots
    mono_models.map(&:uniprot).uniq
  end

  def self.di_uniprots
    di_models.map(&:uniprot).uniq
  end

  def self.all_uniprots
    all_models.map(&:uniprot).uniq
  end

  MonoCollections = ['HL', 'HO', 'SR', 'JA', 'SE', 'SMF', 'SMI', 'CM', 'PAPAM', 'AD']
  DiCollections = ['SDF', 'SDI', 'CD', 'PAPAD', 'DIAD', 'DIHL']
  ChipseqCollections = ['CM', 'CD']
  SelexRebuiltCollections = ['SMF', 'SMI', 'SDF', 'SDI']


  def self.hocomoco_qualities
    @hocomoco_qualities ||= File.readlines('hocomoco_qualities.tsv').map{|line|
      line.chomp.split("\t")
    }.to_h
  end

  # Models to take when no validation supplied in order of reliability
  MonoCollectionsReliability = [
    'PAPAM',
    ->(model){
      model.collection_short_name == 'HL' && ['A', 'B', 'C'].include?( hocomoco_qualities[model.model_name] )
    },
    'CM',
    ->(model){
      model.collection_short_name == 'HL' && hocomoco_qualities[model.model_name] == 'D'
    },
    'SMI', 'SMF']
  DiCollectionsReliability = ['PAPAD']

  # We take only denovo collections and hocomoco legacy for a final bundle
  CollectionsForFinalBundle = ChipseqCollections + SelexRebuiltCollections + ['HL', 'PAPAM', 'PAPAD']

  MonoCollectionsForFinalBundle = CollectionsForFinalBundle & MonoCollections
  DiCollectionsForFinalBundle = CollectionsForFinalBundle & DiCollections

  ## We take only integrated collections from SELEX rebuilt when validation not available
  ## And we don't take any SELEX motifs for dinucleotide models
  ##
  ## CollectionsForFinalBundleWithoutValidation - has an actual but unused value
  ## (because ordered lists are used instead of sets to prioritize models by collection)
  # CollectionsForFinalBundleWithoutValidation = ['HL', 'CM', 'CD', 'SMI']

  CollectionNames = {
    'HL' => 'HOCOMOCO-v10',
    'DIHL' => 'DI-HOCOMOCO-v10',
    'HO' => 'HOMER',
    'JA' => 'JASPAR',
    'SR' => 'SWISSREGULON',
    'SE' => 'HTSELEX',
    'AD' => 'ADDITIONAL-EXTERNAL-MOTIFS',
    'DIAD' => 'ADDITIONAL-EXTERNAL-DI-MOTIFS',

    'CM' => 'CHIPSEQ',
    'SMF' => 'HTSELEX-R',
    'SMI' => 'HTSELEX-I',

    'CD' => 'DI-CHIPSEQ',
    'SDF' => 'DI-HTSELEX-R',
    'SDI' => 'DI-HTSELEX-I',

    'PAPAM' => 'PLURIP',
    'PAPAD' => 'DI-PLURIP',
  }
end
