require_relative 'symmetric_matrix'
require_relative 'distance_matrix'

module Clusterize
  class SimilarityMatrix < SymmetricMatrix
    def initialize(lower_triangle, items)
      super
      raise Error, 'Similarities should be in range [0;1]'  unless lower_triangle.all?{|row| row.all?{|el| (0.0 .. 1.0).include?(el) } }
      raise Error, 'Item should be similar (similarity = 1.0) to itself'  unless size.times.all?{|index| lower_triangle[index][index] == 1 }
    end

    def to_distance_matrix
      lower_triangle_distances = lower_triangle.map{|row| row.map{|el| 1.0 - el } }
      DistanceMatrix.new(lower_triangle_distances, items)
    end

    def ==(other)
      other.is_a?(SimilarityMatrix) && super
    end
  end
end
