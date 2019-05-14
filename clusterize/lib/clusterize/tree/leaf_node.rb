require 'set'
require_relative 'tree_node'
require_relative '../errors'

module Clusterize
  module Tree
    class LeafNode

      include TreeNode
      attr_reader :parent, :data

      def initialize(parent, data = {})
        @parent, @data = parent, data
      end

      def children
        Set.new
      end

      def leaf?
        true
      end

      def to_s(&block)
        info = block_given? ? block.call(self) : data.to_s
        "`#{info}`"
      end

      alias_method :inspect, :to_s

      def each_node(&block)
        return enum_for(:each_node)  unless block_given?
        block.call(self)
      end

      def each_leaf(&block)
        return enum_for(:each_leaf)  unless block_given?
        block.call(self)
      end

      # def each_node_with_parent(options = {}, &block)
      #   return enum_for(:each_node_with_parent, options)  unless block_given?
      #   block.call(self, options[:parent])
      # end

      def number_of_leafs
        1
      end

      def number_of_nodes
        1
      end

      def subtree_nodes
        [self]
      end

      def subtree_leafs
        [self]
      end

      # def map_tree(&block)
      #   Leaf.new(block.call(self))
      # end

      def ==(other)
        other.is_a?(LeafNode) && data == other.data && parent == parent
      end

      alias_method :eql?, :==

      def hash
        [data, parent].hash
      end

      # def subtree_clusters(&split_children_block)
      #   return [[self]]
      # end

      # def subtree_clusters(&split_children_block)
      #   return Forest.new([self])
      # end

    end
  end
end
