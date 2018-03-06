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

  def auc_infos(all_aucs)
    result = models.select{|model| all_aucs.has_key?(model) }.map{|model|
      [model, all_aucs[model]]
    }.to_h
    AUCs.new(result)
  end
end
