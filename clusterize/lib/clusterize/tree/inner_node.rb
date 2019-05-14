# TODO:
# 1) Isn't #hash method too slow because of recursion. What to do with it? Caching? Looking only for nearest children (what to do with no-data nodes in such case?)?
# 2) It's misleading that Leaf isn't Tree. What about other names, common ancestor etc?
# 3) splitting into subtrees

require 'set'
require_relative 'tree_node'
require_relative '../errors'

module Clusterize
  module Tree
    class InnerNode

      include TreeNode
      attr_reader :parent, :data, :children

      def initialize(parent, children, data = {})
        @parent, @data = parent, data
        children_set = Set.new(children).freeze
        raise Error, 'One can\'t create Tree node without children, use Leaf'  if children_set.empty?
        raise Error, 'Node children should be unique'  unless children.size == children_set.size
        @children = children_set
      end

      def leaf?
        false
      end

      def to_s(&block)
        info = block_given? ? block.call(self) : data.to_s
        "(" + children.map{|child| child.to_s(&block) }.join(',') + "):`#{info}`"
      end

      alias_method :inspect, :to_s

      # Depth First Search order
      def each_node(&block)
        return enum_for(:each_node)  unless block_given?
        @children.each do |child|
          child.each_node(&block)
        end
        block.call(self)
      end

      def each_leaf(&block)
        return enum_for(:each_leaf)  unless block_given?
        @children.each do |child|
          child.each_leaf(&block)
        end
      end

      def number_of_nodes
        @number_of_nodes ||= @children.map(&:number_of_nodes).inject(1, &:+)
      end

      def number_of_leafs
        @number_of_leafs ||= @children.map(&:number_of_leafs).inject(0, &:+)
      end
      
      def subtree_nodes
        @children.map(&:subtree_nodes).inject([self], &:+)
      end

      def subtree_leafs
        @children.map(&:subtree_leafs).inject([], &:+)
      end

      # def map_tree(&block)
      #   self.class.new(children.map{|node| node.map_tree(&block) }, block.call(self))
      # end

      def ==(other)
        other.is_a?(InnerNode) && data == other.data && children == other.children && parent == parent
      end

      alias_method :eql?, :==

      def hash
        [data, parent, children].hash
      end

  ## Is it necessary?

      # # ATTENTION! This method's block now has inverse meaning
      # # glues child nodes until we need to set a split point. All upward nodes of splitted are in separate clusters
      # # (((1,2),(3,4)),((5,6),(7,8))) with block ->(node){node.children.max >=8 } -->  [[1,2,3,4],[5,6],[7],[8]]
      # def subtree_clusters(&split_children_block)
      #   # split_children_block
      #   raise '#subtree_clusters needs a block to select splitting points'  unless block_given?
      #   clusters_for_branches = children.map{|child| child.subtree_clusters(&split_children_block) }
      #   if split_children_block.call(self) || clusters_for_branches.any?{|branch_clusters| branch_clusters.size > 1 }
      #     branch_clusters.inject(&:+)  # [[1],[2]] + [[3]] --> [[1],[2].[3]]
      #   else
      #     [clusters_for_branches.map(&:first).inject(&:+)]  # [[1,2]] + [[3]] --> [[1,2,3]]
      #   end
      # end

      # def subtree_clusters(&split_children_block)
      #   # split_children_block
      #   raise '#subtree_clusters needs a block to select splitting points'  unless block_given?
      #   child_forests = children.map{|child| child.subtree_clusters(&split_children_block) }
      #   if split_children_block.call(self) || child_forests.any?{|branch_clusters| branch_clusters.number_of_trees > 1 }
      #     child_forests.inject(Forest.new, &:+)  # [[1],[2]] + [[3]] --> [[1],[2].[3]]
      #   else
      #     Forest.new([self]) # [[1,2]] + [[3]] --> [[1,2,3]]
      #   end
      # end
    end
  end
end
