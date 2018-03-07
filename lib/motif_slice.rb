class MotifsSlice
  attr_reader :slice_type, :semiuniprot
  attr_reader :chipseq_motifs, :hocomoco_motifs, :additional_motifs
  def initialize(slice_type:, semiuniprot:,
                chipseq_motifs:, hocomoco_motifs:, additional_motifs:)
    @slice_type = slice_type
    @semiuniprot = semiuniprot
    @chipseq_motifs = chipseq_motifs
    @hocomoco_motifs = hocomoco_motifs
    @additional_motifs = additional_motifs
  end

  def species_with_currated_motifs
    chipseq_motifs.map{|motif|
      motif.split('.').first.split('_').last
    }.uniq.sort
  end

  def self.from_file(slice_fn, model_type)
    case model_type.to_sym
    when :mono
      MonoMotifsSlice.from_file(slice_fn)
    when :di
      DiMotifsSlice.from_file(slice_fn)
    else
      raise "Undefined model type #{model_type}"
    end
  end
end
