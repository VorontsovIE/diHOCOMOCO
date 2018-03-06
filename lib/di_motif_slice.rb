require_relative 'motif_slice'

class DiMotifsSlice < MotifsSlice
  def self.banlist
    File.readlines('curation/banlist_di.txt').map(&:strip)
  end

  def self.all_hocomoco_motifs(semiuniprot)
    Dir.glob("models/pcm/di/hocomoco_legacy/#{semiuniprot}_*.H10DI.*.dpcm").map{|fn|
      File.basename(fn, '.dpcm')
    }.reject{|motif|
      self.banlist.include?(motif)
    }
  end

  def self.from_file(slice_fn)
    motif_final = File.basename(slice_fn, '.txt')
    slice_type = motif_final.split('.').last[0]
    semiuniprot = motif_final.split('.').first  # without species part
    all_motifs = File.readlines(slice_fn).map(&:chomp)

    additional_motifs = all_motifs.select{|motif|
      motif.match(/~DIAD~/)
    }.map{|motif| motif.split('~').last }

    motifs = all_motifs.reject{|motif|
      motif.match(/~DIAD~/)
    }
    hocomoco_motifs = motifs.select{|motif|
      motif.match(/\.H10DI\./)
    }
    chipseq_motifs = motifs.reject{|motif|
      motif.match(/\.H10DI\./)
    }.reject{|motif|
      banlist.include?(motif)
    }

    hocomoco_motifs = self.all_hocomoco_motifs(semiuniprot)  if hocomoco_motifs.empty? && slice_type == 'M'

    self.new( slice_type: slice_type, semiuniprot: semiuniprot,
              chipseq_motifs: chipseq_motifs,
              hocomoco_motifs: hocomoco_motifs,
              additional_motifs: additional_motifs )
  end

  def models
    @models_cache ||= begin
      result = []
      result += chipseq_motifs.map{|motif_name|
        uniprot = motif_name.split('.').first
        Model.new("#{uniprot}~CD~#{motif_name}", :di)
      }
      result += hocomoco_motifs.map{|motif_name|
        uniprot = motif_name.split('.').first
        Model.new("#{uniprot}~DIHL~#{motif_name}", :di)
      }
      result += additional_motifs.map{|motif_name|
        uniprot = motif_name.split('.').first
        Model.new("#{uniprot}~DIAD~#{motif_name}", :di)
      }
      result
    end
  end
end
