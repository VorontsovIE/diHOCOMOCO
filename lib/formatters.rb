require 'bioinform'
require 'information_content'

def pcm2ppm(pcm)
  pcm.map{|pos|
    count = pos.inject(0.0, &:+)
    pos.map {|el| el / count }
  }
end

def in_meme_format(infos_for_motifs)
  result = ""
  result << "MEME version 4\n\n"
  result << "ALPHABET= ACGT\n\n"
  result << "strands: + -\n\n"
  result << "Background letter frequencies\n"
  result << "A 0.25 C 0.25 G 0.25 T 0.25\n\n"

  infos_for_motifs.each do |motif_infos|
    count = motif_infos[:pcm][0].inject(&:+).round
    ppm = pcm2ppm(motif_infos[:pcm])
    result << "MOTIF #{motif_infos[:name]}\n"
    result << "letter-probability matrix: alength= 4 w= #{motif_infos[:length]} nsites= #{count}\n"
    result << ppm.map{|pos| pos.join("\t") }.join("\n") << "\n"
    result << "URL http://hocomoco.autosome.ru/motif/#{motif_infos[:name]}\n"
    result << "\n"
  end
  result
end

def in_transfac_format(infos_for_motifs)
  result = ""
  result << "VV  HOCOMOCOv11, SEP 2017\n"
  result << "XX\n"
  result << "//\n"
  infos_for_motifs.each_with_index do |motif_infos, index|
    consensus = motif_infos[:consensus_string]
    rounded_pcm = round_pcm(motif_infos[:pcm])

    uniprot = motif_infos[:uniprot]
    species = motif_infos[:species].downcase

    result << "AC  M" << ('%05d' % (index + 1001)) << "\n" # offset 1000 is taken to distinguish HOCOMOCO v9 from v10
    result << "XX\n"
    result << "ID  #{motif_infos[:name]}\n"
    result << "XX\n"
    result << "NA  #{motif_infos[:uniprot]}\n"
    result << "XX\n"
    if species == 'human'
      result << "BF  #{motif_infos[:uniprot]}; Species: human, Homo sapiens\n"
    elsif species == 'mouse'
      result << "BF  #{motif_infos[:uniprot]}; Species: mouse, Mus musculus\n"
    else
      result << "BF  #{motif_infos[:uniprot]}; Species: #{species}\n"
    end
    result << "XX\n"
    result << "P0      A      C      G      T\n"
    rounded_pcm.each_with_index do |pos, pos_index|
      result << "%02d %6.20g %6.20g %6.20g %6.20g      %s\n" % [pos_index + 1, *pos, consensus[pos_index]]
    end
    result << "XX\n"
    result << "//\n"
  end
  result
end

def rounded_pcm_position(pos, expected_count)
  rounded_pos = pos.map(&:round)
  diff = expected_count - rounded_pos.inject(&:+)
  order = pos.each_index.sort{|ind_1, ind_2|
    rounded_pos[ind_1] <=> rounded_pos[ind_2]
  }
  if diff > 0
    rounded_pos[order.first] += diff # smallest count increase
  elsif diff < 0
    rounded_pos[order.last] += diff # largest count decreased
  end

  raise 'Unable to round matrix'  unless rounded_pos.inject(&:+) == expected_count
  rounded_pos
end

def round_pcm(pcm_matrix)
  counts = pcm_matrix.map{|pos| pos.inject(&:+) }
  raise 'Different counts'  unless counts.all?{|count| (count - counts[0]).to_f / counts[0]  < 0.001 }
  rounded_word_count = counts[0].round

  pcm_matrix.map{|pos|
    rounded_pcm_position(pos, rounded_word_count)
  }
end

def in_jaspar_format(infos_for_motifs)
  infos_for_motifs.map{|motif_infos|
    pcm_rounded = round_pcm(motif_infos[:pcm])
    matrix_str = pcm_rounded.transpose.map{|nucleotide_counts| nucleotide_counts.join("\t") }.join("\n")
    ">#{motif_infos[:name]}\n#{matrix_str}"
  }.join("\n")
end

def in_homer_format(infos_for_motifs, pvalue:)
  infos_for_motifs.map{|motif_infos|

    #  HOMER always works with uniform background substitution
    #(see "Motif Scanning" section at http://homer.salk.edu/homer/motif/creatingCustomMotifs.html)
    header = [">#{motif_infos[:consensus_string]}", motif_infos[:name], motif_infos[:standard_thresholds][pvalue]].join("\t")

    pseudocount_calc = ->(pos){
      max_count = pos.inject(0.0, &:+)
      Math.log([max_count, 2].max)
    }
    # necessary only for standard pseudocount calculation
    logodds_converter = Bioinform::ConversionAlgorithms::PCM2PWMConverterDifferentCount.new(pseudocount: pseudocount_calc)

    ppm_matrix_corrected = motif_infos[:pcm].map{|pos|
      count = pos.inject(0.0, &:+)
      pseudocount = logodds_converter.calculate_pseudocount(pos)
      pos.map{|el| (el + pseudocount * 0.25) / (count + pseudocount) }
    }

    matrix_str = ppm_matrix_corrected.map{|pos| pos.join("\t") }
    [header, *matrix_str].join("\n")
  }.join("\n")
end
