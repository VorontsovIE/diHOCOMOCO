require 'bioinform'

##########
module Bioinform
  module MotifModel
    class DiPM # Doesn't work with alphabet
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
            Math.log(numerator / denominator)
          }
        }

        DiPWM.new(converted_matrix)
      end
    end

    class DiPWM < DiPM
    end
  end
end


# str = File.read('models/pcm/di/all/ALX1_HUMAN~SDF~ALX1_HUMAN_DBD_0.dpcm')
# parser = Bioinform::MatrixParser.new(fix_nucleotides_number: 16)
# infos = parser.parse(str)
# pcm = Bioinform::MotifModel::DiPCM.new(infos[:matrix]).named(infos[:name])
# pwm = pcm.to_pwm
# puts pwm
