class FastaSequence
  attr_reader :sequence, :header
  def initialize(sequence, header)
    @sequence = sequence.upcase
    @header = header
  end

  def length
    sequence.length
  end

  # yields nucleotide and weight for each position
  def each_position(&block)
    sequence.each_char(&block)
  end

  def self.each_in_file(filename)
    return enum_for(:each_in_file, filename)  unless block_given?
    File.open(filename) do |f|
      f.each_line.map(&:chomp).each_slice(2) do |header, sequence|
        yield self.new(sequence, header)
      end
    end
  end
end
