class WeightedSequence
  attr_reader :sequence, :weights
  def initialize(sequence, weights)
    raise 'Weights are incompatible with sequence'  unless weights.length == sequence.length
    @sequence = sequence
    @weights = weights
  end

  def length
    sequence.length
  end

  # yields nucleotide and weight for each position
  def each_position(&block)
    sequence.each_char.zip(weights).each(&block)
  end

  def self.each_in_file(filename)
    return enum_for(:each_in_file, filename)  unless block_given?
    File.open(filename) do |f|
      f.each_line.map(&:chomp).each_slice(2) do |weights_line, sequence|
        weights = weights_line.sub(/^>\s*/, '').split.map(&:to_f)
        yield self.new(sequence, weights)
      end
    end
  end
end
