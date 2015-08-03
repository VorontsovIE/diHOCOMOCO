require_relative 'weighted_sequence'
require_relative 'counts'
require_relative 'background'

class SequenceDataset
  attr_reader :filename

  # control/control/AEBP2_HUMAN^PEAKS030225.mfa
  def initialize(filename)
    @filename = filename
  end

  def name
    filename.pathmap('%n')
  end

  def uniprot
    name.split('^').first
  end

  def each_sequence
    return enum_for(:each_sequence)  unless block_given?
    WeightedSequence.each_in_file(filename){|weighted_sequence|
      yield weighted_sequence
    }
  end

  def local_mono_background
    @local_mono_background ||= begin
      counts = Hash.new(0)
      each_sequence{|weighted_sequence|
        weighted_sequence.each_position{|letter, weight|
          counts[letter] += 1 # weight
        }
      }
      MonoCounts.from_hash(counts).plus_revcomp.frequencies
    end
  end

  def local_di_background
    @local_di_background ||= begin
      counts = Hash.new(0)
      each_sequence{|weighted_sequence|
        weighted_sequence.each_position.each_cons(2){|(letter_1, weight_1), (letter_2, weight_2)|
          diletter = "#{letter_1}#{letter_2}"
          counts[diletter] += 1 # (weight_1 + weight_2) / 2.0
        }
      }
      DiCounts.from_hash(counts).plus_revcomp.frequencies
    end
  end

  def self.each_by_glob(glob)
    return enum_for(:each_by_glob, glob)  unless block_given?
    FileList[glob].each{|filename|
      yield SequenceDataset.new(filename)
    }
  end
end
