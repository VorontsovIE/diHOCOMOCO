require_relative 'fasta_sequence'
require_relative 'counts'
require_relative 'frequencies'

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
    name.split('.').first
  end

  def mono_models; Models.mono_models_by_uniprot(uniprot); end
  def di_models; Models.di_models_by_uniprot(uniprot); end
  def all_models; Models.all_models_by_uniprot(uniprot); end

  def each_sequence
    return enum_for(:each_sequence)  unless block_given?
    FastaSequence.each_in_file(filename){|fasta_sequence|
      yield fasta_sequence
    }
  end

  # def local_mono_background_path
  #   "control/local_backgrounds/mono/#{name}.txt"
  # end

  def local_di_background_path
    "control/local_backgrounds/di/#{name}.txt"
  end

  # def local_mono_background
  #   @local_mono_background ||= begin
  #     counts = Hash.new(0)
  #     each_sequence{|fasta_sequence|
  #       fasta_sequence.each_position{|letter|
  #         counts[letter] += 1
  #       }
  #     }
  #     MonoCounts.from_hash(counts).plus_revcomp.frequencies
  #   end
  # end

  def local_di_background
    @local_di_background ||= begin
      counts = Hash.new(0)
      each_sequence{|fasta_sequence|
        fasta_sequence.each_position.each_cons(2){|letter_1, letter_2|
          diletter = "#{letter_1}#{letter_2}"
          counts[diletter] += 1
        }
      }
      DiCounts.from_hash(counts).plus_revcomp.frequencies
    end
  end

  def self.each_file_by_glob(glob)
    return enum_for(:each_file_by_glob, glob)  unless block_given?
    FileList[glob].each{|filename|
      yield self.new(filename)
    }
  end

  def self.each_dataset(&block)
    each_file_by_glob('control/control/*.mfa', &block)
  end

  def self.each_for_uniprot(uniprot, &block)
    each_file_by_glob("control/control/#{uniprot}.*.mfa", &block)
  end

  def self.each_uniprot(&block)
    each_dataset.map(&:uniprot).uniq.sort.each(&block)
  end

  def self.by_name(dataset_name)
    self.new("control/control/#{dataset_name}.mfa")
  end
end
