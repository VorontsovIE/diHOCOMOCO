require 'bioinform'

##########
class Bioinform::MotifModel::PCM
  private def rounded_position(pos, expected_count)
    rounded_pos = pos.map(&:round)
    diff = expected_count - rounded_pos.inject(&:+)
    order = pos.each_index.sort{|ind_1, ind_2|
      rounded_pos[ind_1] <=> rounded_pos[ind_2]
    }
    if diff > 0
      rounded_pos[order.first] += diff # smallest count increase
    elsif diff < 0
      rounded_pos[order.last] += diff # largest count decreased
    end

    raise 'Unable to round matrix'  unless rounded_pos.inject(&:+) == expected_count
    rounded_pos
  end

  def round
    counts = each_position.map{|pos| pos.inject(&:+) }
    raise 'Different counts'  unless counts.all?{|count| (count - counts[0]).to_f / counts[0]  < 0.001 }
    rounded_word_count = counts[0].round

    rounded_matrix = each_position.map{|pos|
      rounded_position(pos, rounded_word_count)
    }
    self.class.new(rounded_matrix)
  end
end
