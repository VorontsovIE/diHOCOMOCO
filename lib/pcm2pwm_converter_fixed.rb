require 'bioinform'

module Bioinform
  module ConversionAlgorithms
    class PCM2PWMConverterDifferentCount

      attr_reader :background, :pseudocount
      def initialize(options = {})
        @background = options.fetch(:background, Bioinform::Background::Uniform)
        @pseudocount = options.fetch(:pseudocount, :log)
      end

      def calculate_pseudocount(pos)
        count = pos.inject(0.0, &:+)
        case @pseudocount
        when Numeric
          @pseudocount
        when :log
          Math.log([count, 2].max)
        when :sqrt
          Math.sqrt(count)
        when Proc
          @pseudocount.call(pos)
        else
          raise Error, 'Unknown pseudocount type use numeric or :log or :sqrt or Proc with taking pcm parameter'
        end
      end


      def convert(pcm)
        raise Error, "#{self.class}#convert accepts only models acting as PCM"  unless MotifModel.acts_as_pcm?(pcm)
        conv_matrix = pcm.each_position.map do |pos|
          actual_pseudocount = calculate_pseudocount(pos)
          count = pos.inject(0.0, &:+)
          pos.each_index.map do |letter|
            numerator = pos[letter] + @background.frequencies[letter] * actual_pseudocount
            denominator = @background.frequencies[letter]*(count + actual_pseudocount)
            Math.log(numerator.to_f / denominator)
          end
        end

        pwm = MotifModel::PWM.new(conv_matrix)
        if pcm.respond_to? :name
          pwm.named(pcm.name)
        else
          pwm
        end
      end
    end
  end
end


# хорошо бы в macro-perfectos-ape JAVA встроить
def pcm_to_pwm(from_file, to_file)
  parser = Bioinform::MatrixParser.new(fix_nucleotides_number: 4)
  infos = parser.parse(File.read(from_file))
  name = infos[:name] || File.basename(from_file, '.pcm')
  pcm = Bioinform::MotifModel::PCM.new(infos[:matrix]).named(name)

  pseudocount_calc = ->(pos){
    max_count = pos.inject(0.0, &:+)
    Math.log([max_count, 2].max)
  }
  converter = Bioinform::ConversionAlgorithms::PCM2PWMConverterDifferentCount.new(pseudocount: pseudocount_calc)
  pwm = converter.convert(pcm)

  File.write to_file, pwm.to_s
end
