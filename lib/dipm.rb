require 'bioinform'
require 'fileutils'

##########
module Bioinform
  module MotifModel
    class DiPM # Doesn't work with alphabet

      Nucleotides = ['A', 'C', 'G', 'T']
      Dinucleotides = Nucleotides.product(Nucleotides).map(&:join)

      RevcompIndexMapping = Dinucleotides.each_with_index.map{|dinucleotide, dinucleotide_index|
        reverse_dinucleotide = dinucleotide.reverse.tr('ACGT','TGCA')
        [dinucleotide_index, Dinucleotides.index(reverse_dinucleotide)]
      }.to_h

      attr_reader :matrix
      def initialize(matrix)
        @matrix = matrix
        raise ValidationError.new('invalid matrix', validation_errors: validation_errors)  unless valid?
      end

      def validation_errors
        errors = []
        errors << "matrix should be an Array"  unless matrix.is_a? Array
        errors << "matrix shouldn't be empty"  unless matrix.size > 0
        errors << "each matrix position should be an Array"  unless matrix.all?{|pos| pos.is_a?(Array) }
        errors << "each matrix position should be of size compatible with alphabet (=#{16})"  unless matrix.all?{|pos| pos.size == 16 }
        errors << "each matrix element should be Numeric"  unless matrix.all?{|pos| pos.all?{|el| el.is_a?(Numeric) } }
        errors
      end
      private :validation_errors

      def valid?
       validation_errors.empty?
      rescue
        false
      end

      private :valid?

      def to_s
        MotifFormatter.new.format(self)
      end

      def named(name)
        NamedModel.new(self, name)
      end

      def length
        matrix.size + 1
      end

      def ==(other)
        self.class == other.class && matrix == other.matrix # alphabet should be considered (when alphabet implemented)
      end

      def each_position
        if block_given?
          matrix.each{|pos| yield pos}
        else
          self.to_enum(:each_position)
        end
      end

      def revcomp
        revcomp_matrix = matrix.map{|pos|
          pos.each_index.map{|index| pos[RevcompIndexMapping[index]] }
        }.reverse
        self.class.new(revcomp_matrix)
      end

      def to_hash
        Dinucleotides.zip(matrix.transpose).to_h
      end

      def self.sum_by_second_letter(pos)
        4.times.map{|first_letter|
          4.times.map{|second_letter|
            pos[4*first_letter + second_letter]
          }.inject(0.0, &:+)
        }
      end

      def self.sum_by_first_letter(pos)
        4.times.map{|second_letter|
          4.times.map{|first_letter|
            pos[4*first_letter + second_letter]
          }.inject(0.0, &:+)
        }
      end
    end

    class DiPCM < DiPM
      # canonicaly this is done via converter, but it's too complicated to be implemented right now
      # dibackground also should be a class and so on... Sometime later I'll refactor it
      def to_pwm(background = Array.new(16, 0.0625))
        converted_matrix = @matrix.map{|pos|
          count = pos.inject(0, &:+)
          pseudocount = Math.log([count, 2].max)
          pos.each_index.map{|letter|
            numerator = pos[letter] + background[letter] * pseudocount
            denominator = background[letter] * (count + pseudocount)
            Math.log(numerator.to_f / denominator)
          }
        }

        DiPWM.new(converted_matrix)
      end

      def to_mono
        matrix_prefix = matrix.map{|pos| DiPM.sum_by_second_letter(pos) }
        matrix_last_pos = DiPM.sum_by_first_letter(matrix.last)
        PCM.new(matrix_prefix + [matrix_last_pos])
      end
    end

    class DiPWM < DiPM
    end
  end
end


# хорошо бы в macro-perfectos-ape JAVA встроить
def dipcm_to_dipwm(from_file, to_file)
  output_dir = File.dirname(to_file)
  FileUtils.mkdir_p(output_dir)  unless Dir.exist?(output_dir)
  parser = Bioinform::MatrixParser.new(fix_nucleotides_number: 16, nucleotides_in: :columns)
  infos = parser.parse(File.read(from_file))
  name = infos[:name] || File.basename(from_file, '.dpcm')
  pcm = Bioinform::MotifModel::DiPCM.new(infos[:matrix]).named(name)
  pwm = pcm.to_pwm
  File.write to_file, pwm.to_s
end
