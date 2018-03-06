require_relative 'motif_slice'

class MonoMotifsSlice < MotifsSlice
  def self.banlist
    File.readlines('curation/banlist_mono.txt').map(&:strip)
  end

  def self.all_hocomoco_motifs(semiuniprot)
    Dir.glob("models/pcm/mono/hocomoco_legacy/#{semiuniprot}_*.H10MO.*.pcm").map{|fn|
      File.basename(fn, '.pcm')
    }.reject{|motif|
      banlist.include?(motif)
    }
  end

  def self.from_file(slice_fn)
    motif_final = File.basename(slice_fn, '.txt')
    slice_type = motif_final.split('.').last[0]
    semiuniprot = motif_final.split('.').first  # without species part
    all_motifs = File.readlines(slice_fn).map(&:chomp)

    additional_motifs = all_motifs.select{|motif|
      motif.match(/~AD~/)
    }.map{|motif| motif.split('~').last }

    motifs = all_motifs.reject{|motif|
      motif.match(/~AD~/)
    }
    hocomoco_motifs = motifs.select{|motif|
      motif.match(/\.H10MO\./)
    }
    chipseq_motifs = motifs.reject{|motif|
      motif.match(/\.H10MO\./)
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
        Model.new("#{uniprot}~CM~#{motif_name}", :mono)
      }
      result += hocomoco_motifs.map{|motif_name|
        uniprot = motif_name.split('.').first
        Model.new("#{uniprot}~HL~#{motif_name}", :mono)
      }
      result += additional_motifs.map{|motif_name|
        uniprot = motif_name.split('.').first
        Model.new("#{uniprot}~AD~#{motif_name}", :mono)
      }
      result
    end
  end
end
