require_relative 'symmetric_matrix'

# TODO: triangle inequality check
module Clusterize
  class DistanceMatrix < SymmetricMatrix
    def initialize(lower_triangle, items)
      super
      raise Error, 'Distances should be non-negative'  unless lower_triangle.all?{|row| row.all?{|el| el >= 0.0 } }
      raise Error, 'Item have non-null distance to itself'  unless size.times.all?{|index| lower_triangle[index][index] == 0 }
    end

    def to_distance_matrix
      self
    end

    # find index (i,j) corresponding to a minimal distance between items
    # return value ci,cj is always ordered: ci < cj
    def index_of_minimal_element(indices_subset = size.times.to_a)
      raise Error, 'Distance matrix has the only item so looking for two nearest items not possible'  if indices_subset.size < 2
      each_pair(indices_subset).min_by{|ci,cj| element_at(ci, cj) }.sort.reverse
    end

    def each_pair(indices_subset = size.times, &block)
      indices_subset.combination(2, &block)
    end

    # distance matrix with one more node, whose distances are given as `node_distances`
    def with_node(node_distances, node_item)
      raise Error, 'Distances to all items should be of the same size as number of items'  unless node_distances.size == size
      node_distances_expanded = node_distances + [0.0] # distance of node to itself included
      lower_triangle_expanded = lower_triangle + [node_distances_expanded]
      items_expanded = items + [node_item]
      DistanceMatrix.new(lower_triangle_expanded, items_expanded)
    end

    def ==(other)
      other.is_a?(DistanceMatrix) && super
    end
  end
end
