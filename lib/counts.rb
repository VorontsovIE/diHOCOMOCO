require_relative 'support'
require_relative 'frequencies'

# abstract class
class CountsBase
  def self.possible_letters
    raise NotImplemented, 'Should be implemented in a subclass'
  end

  def self.create_frequencies(frequencies)
    raise NotImplemented, 'Should be implemented in a subclass'
  end

  def create_frequencies(frequencies)
    self.class.create_frequencies(frequencies)
  end

  def possible_letters
    self.class.possible_letters
  end

  attr_reader :counts
  def initialize(counts)
    raise "Counts should be an Array with #{possible_letters.size}"  unless counts.is_a?(Array) && counts.size == possible_letters.size
    @counts = counts
  end

  def each_count_with_letter(&block)
    @counts.zip(possible_letters).each(&block)
  end

  def total_count
    @total_count ||= @counts.inject(0.0, &:+)
  end

  def frequencies
    create_frequencies( @counts.map{|count| count.to_f / total_count } )
  end

  private def revcomp(seq)
    seq.tr('ACGTacgt','TGCAtgca').reverse
  end

  # adds counts for reverse complement keys
  def plus_revcomp
    counts_plus_revcomp = Hash.new(0)
    each_count_with_letter{|count, letter| # letter can be diletter
      counts_plus_revcomp[letter] += count
      counts_plus_revcomp[revcomp(letter)] += count
    }
    self.class.from_hash(counts_plus_revcomp)
  end

  def self.from_string(str)
    self.new(str.strip.split(/\s+|,/).map(&:to_f))
  end

  extend InitializeFromHash
  # ignores N-s and other non-ACGT keys
  def self.from_hash(hsh)
    self.new(normalized_hash(hsh).values_at(*self.possible_letters))
  end
end

def Counts(letters, frequencies_base)
  Class.new(CountsBase) do
    singleton_class.send(:define_method, :possible_letters){ letters }
    singleton_class.send(:define_method, :create_frequencies){|frequencies| frequencies_base.new(frequencies) }
  end
end

MonoCounts = Counts(Nucleotides, MonoFrequencies)
DiCounts = Counts(Dinucleotides, DiFrequencies)
