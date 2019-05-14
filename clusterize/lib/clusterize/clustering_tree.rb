require_relative 'matrix/distance_matrix'
require_relative 'tree'
require_relative 'errors'

module Clusterize

  class ClusteringTree
    class Leaf < Clusterize::Leaf
      def to_s
        "`#{data[:item]}`"
      end
    end

    class Tree < Clusterize::Tree
      def to_s
        '(' + children.map(&:to_s).join(',') + ')'
      end
    end

    attr_reader :distance_matrix, :clustering_forest
    def initialize(distance_matrix, clustering_forest)
      raise Error, 'Distance matrix and clustering forest size are not equal'  unless distance_matrix.size == clustering_forest.number_of_nodes
      @distance_matrix, @clustering_forest = distance_matrix, clustering_forest
    end

    def self.from_distance_matrix(distance_matrix)
      nodes = distance_matrix.items.map.with_index{|item, index| Leaf.new(item: item, index: index) }
      forest = Forest.new(nodes)
      raise Error, 'Distance matrix has duplicated items'  unless forest.number_of_nodes == nodes.size
      self.new(DistanceMatrix.new(distance_matrix.matrix, nodes), forest)
    end

    def items
      distance_matrix.items
    end

    def top_level_items
      clustering_forest.trees
    end

    def index_of_minimal_element
      distance_matrix.index_of_minimal_element( top_level_items.map{|root_node| root_node.data[:index]} )
    end

    def distances_to_cluster(first_node_index, second_node_index, &distance_to_cluster_union)
      distance_matrix.size.times.map{|i| distance_to_cluster_union.call(distance_matrix, i, first_node_index, second_node_index) }
    end

    def with_link(first_node_index, second_node_index, &distance_to_cluster_union)
      node_distances = distances_to_cluster(first_node_index, second_node_index, &distance_to_cluster_union)
      first_node = items[first_node_index]
      second_node = items[second_node_index]
      data = {item: "Union (#{first_node_index},#{second_node_index})", index: distance_matrix.size}
      union_node = Tree.new([first_node, second_node], data)
      new_distance_matrix = distance_matrix.with_node(node_distances, union_node)
      new_clustering_forest = clustering_forest.with_link([first_node, second_node], join_node: union_node)
      ClusteringTree.new(new_distance_matrix, new_clustering_forest)
    end

    def clusterize(&distance_to_cluster_union)
      current_clustering_tree = self
      while current_clustering_tree.clustering_forest.number_of_trees > 1
        nodes_to_unite = current_clustering_tree.index_of_minimal_element
        current_clustering_tree = current_clustering_tree.with_link(*nodes_to_unite, &distance_to_cluster_union)
      end
      current_clustering_tree
    end

    def cluster_tree(&distance_to_cluster_union)
      clusterize(&distance_to_cluster_union).clustering_forest.trees.first
    end

    def cluster_tree_with_distances(&distance_to_cluster_union)
      cluster_tree(&distance_to_cluster_union).map_tree { |node|
        node.data.merge(distance: nil)
      }
    end
  end
end
