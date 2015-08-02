require_relative 'support'

# abstract class
class FrequenciesBase
  def self.possible_letters
    raise NotImplemented, 'Should be implemented in a subclass'
  end

  def possible_letters
    self.class.possible_letters
  end

  attr_reader :frequencies
  def initialize(frequencies)
    raise "Frequencies should be an Array with #{possible_letters.size} elements"  unless frequencies.is_a?(Array) && frequencies.size == possible_letters.size
    raise 'Sum of frequencies should be 1'  unless (frequencies.inject(0.0, &:+) - 1.0).abs < 1.0e-5
    @frequencies = frequencies
  end

  def to_s
    @frequencies.join(',')
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

def Frequencies(letters)
  Class.new(FrequenciesBase) do
    singleton_class.send(:define_method, :possible_letters){ letters }
  end
end

MonoFrequencies = Frequencies(Nucleotides)
DiFrequencies = Frequencies(Dinucleotides)
