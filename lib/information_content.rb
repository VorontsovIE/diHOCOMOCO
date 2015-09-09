require 'bioinform'

class Float
  def log_fact
    Math.lgamma(self + 1).first
  end
end

class Integer
  def log_fact
    self.to_f.log_fact
  end
end

IUPAC_CODE = {
  'A' => 'A', 'C' => 'C', 'G' => 'G', 'T' => 'T',
  'AG' => 'R', 'CT' => 'Y', 'GT' => 'K', 'AC' => 'M', 'CG' => 'S', 'AT' => 'W',
  'CGT' => 'B', 'AGT' => 'D', 'ACT' => 'H', 'ACG' => 'V',
  'ACGT' => 'N'
}

def infocod(pos)
  words_count = pos.inject(0.0, &:+)
  (pos.map(&:log_fact).inject(0.0, &:+) - words_count.log_fact) / words_count
end

def icd2of4(pos)
  words_count = pos.inject(0.0, &:+)
  i2o4 = words_count / 2.0
  infocod([i2o4, i2o4, 0, 0])
end

def icd3of4(pos)
  words_count = pos.inject(0.0, &:+)
  i3o4 = words_count / 3.0
  infocod([i3o4, i3o4, i3o4, 0])
end
def icdThc(pos); icd3of4(pos); end

def icdTlc(pos)
  words_count = pos.inject(0.0, &:+)
  io = words_count / 6.0
  infocod([2*io, 2*io, io, io])
end

def num_significant_letters(pos)
  icd = infocod(pos)
  if icd > icd2of4(pos)
    1
  elsif icd > icdThc(pos)
    2
  elsif icd > icdTlc(pos)
    3
  else
    4
  end
end

module Bioinform
  module MotifModel
    class PCM
      def consensus_string(beautiful: true)
        nucleotides = ['A', 'C', 'G', 'T']
        each_position.map{|pos|
          scores = pos.uniq.sort.reverse.first(num_significant_letters(pos))

          lets = nucleotides.zip(pos).select{|letter, score|
            scores.include?(score)
          }.map{|letter, score| letter }.inject('', &:+)

          (beautiful && lets.length > 2) ? IUPAC_CODE[lets].downcase : IUPAC_CODE[lets]
        }.inject('', &:+)
      end
    end

    class DiPCM
      def consensus_string(beautiful: true)
        to_mono.consensus_string(beautiful: beautiful)
      end
    end
  end
end
