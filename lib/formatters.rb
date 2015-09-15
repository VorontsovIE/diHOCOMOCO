require 'bioinform'
require 'information_content'
require 'bioinform_support'

def in_meme_format(motif_pcms)
  result = ""
  result << "MEME version 4\n\n"
  result << "ALPHABET= ACGT\n\n"
  result << "strands: + -\n\n"
  result << "Background letter frequencies\n"
  result << "A 0.25 C 0.25 G 0.25 T 0.25\n\n"

  motif_pcms.sort_by(&:name).each do |pcm|
    count = pcm.matrix[0].inject(&:+).round
    ppm = Bioinform::ConversionAlgorithms::PCM2PPMConverter.new.convert(pcm)
    result << "MOTIF #{pcm.name}\n"
    result << "letter-probability matrix: alength= 4 w= #{pcm.length} nsites= #{count}\n"
    result << ppm.matrix.map{|pos| pos.join("\t") }.join("\n") << "\n"
    result << "URL http://hocomoco.autosome.ru/motif/#{pcm.name}\n"
    result << "\n"
  end
  result
end

def in_transfac_format(motif_pcms)
  result = ""
  result << "VV\tHOCOMOCOv10, SEP 2015\n"
  result << "XX\n"
  result << "//\n"
  motif_pcms.sort_by(&:name).each_with_index do |pcm, index|
    consensus = pcm.consensus_string
    rounded_pcm = pcm.round

    uniprot = pcm.name.split('.').first
    species = uniprot.split('_').last.downcase

    result << "AC\tM" << ('%05d' % (index + 1001)) << "\n" # offset 1000 is taken to distinguish HOCOMOCO v9 from v10
    result << "XX\n"
    result << "ID\t#{pcm.name}\n"
    result << "XX\n"
    result << "NA\t#{uniprot}\n"
    result << "XX\n"
    if species == 'human'
      result << "BF\t#{uniprot}; Species: human, Homo sapiens\n"
    elsif species == 'mouse'
      result << "BF\t#{uniprot}; Species: mouse, Mus musculus\n"
    else
      result << "BF\t#{uniprot}; Species: #{species}\n"
    end
    result << "XX\n"
    result << ['P0', 'A', 'C', 'G', 'T'].join("\t") << "\n"
    rounded_pcm.matrix.each_with_index do |pos, pos_index|
      result << ['%02d' % (pos_index + 1), *pos, consensus[pos_index]].join("\t") << "\n"
    end
    result << "XX\n"
    result << "//\n"
  end
  result
end

def in_jaspar_format(motif_pcms)
  motif_pcms.sort_by(&:name).map{|pcm|
    pcm_rounded = pcm.round
    matrix_str = pcm_rounded.matrix.transpose.map{|nucleotide_counts| nucleotide_counts.join("\t") }.join("\n")
    ">#{pcm.name}\n#{matrix_str}"
  }.join("\n")
end

def in_homer_format(motif_pcms, thresholds_by_model, pvalue:)
  motif_pcms.map{|pcm|

    #  псевдокаунты!!!
    header = [">#{pcm.consensus_string}", pcm.name, thresholds_by_model[pcm.name][pvalue]].join("\t")
    matrix_str = pcm.each_position.map{|pos| pos.join("\t") }
    [header, *matrix_str].join("\n")
  }.join("\n")
end
