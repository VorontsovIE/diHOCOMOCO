Nucleotides = ['A', 'C', 'G', 'T'] # We count N's, but don't use them
Dinucleotides = Nucleotides.product(Nucleotides).map(&:join)

class WeightedSequence
  attr_reader :sequence, :weights
  def initialize(sequence, weights)
    raise 'Weights are incompatible with sequence'  unless weights.length == sequence.length
    @sequence = sequence
    @weights = weights
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


def calculate_local_mono_background(filename)
  counts = Hash.new(0)
  WeightedSequence.each_in_file(filename){|weighted_sequence|
    weighted_sequence.each_position{|letter, weight|
      counts[letter] += weight
    }
  }
  # Nucleotides is ACGT; We count N's, but don't use them
  total_count = Nucleotides.map{|letter| counts[letter] }.inject(:+)
  Nucleotides.map{|letter| counts[letter].to_f / total_count }
end

def calculate_local_di_background(filename)
  counts = Hash.new(0)
  WeightedSequence.each_in_file(filename){|weighted_sequence|
    weighted_sequence.each_position.each_cons(2){|(letter_1, weight_1), (letter_2, weight_2)|
      diletter = "#{letter_1}#{letter_2}"
      counts[diletter] += (weight_1 + weight_2) / 2.0
    }
  }
  # Nucleotides is ACGT; We count N's, but don't use them
  total_count = Dinucleotides.map{|diletter| counts[diletter] }.inject(:+)
  Dinucleotides.map{|diletter| counts[diletter].to_f / total_count }
end

desc 'Calculate local background for each dataset'
task :calculate_local_backgrounds do
  mkdir_p 'control/local_backgrounds/mono'
  mkdir_p 'control/local_backgrounds/di'
  FileList['control/control/*.mfa'].each{|fn|
    File.write(fn.pathmap('control/local_backgrounds/mono/%n.txt'), calculate_local_mono_background(fn).join("\t"))
    File.write(fn.pathmap('control/local_backgrounds/di/%n.txt'), calculate_local_di_background(fn).join("\t"))
  }
end
