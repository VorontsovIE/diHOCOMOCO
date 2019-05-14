require 'set'
require_relative 'inner_node'
require_relative 'leaf_node'
require_relative '../errors'

module Clusterize
  class Forest
    attr_reader :root_nodes
    def initialize(root_nodes = [])
      @root_nodes = Set.new(root_nodes)
      raise Error, 'Node children should be unique'  unless @root_nodes.size == root_nodes.size
    end

    def trees
      root_nodes
    end

    def to_s
      '{' + root_nodes.map(&:to_s).join('; ') + '}'
    end
    alias_method :inspect, :to_s

    def with_link(nodes_to_unite, options = {})
      raise Error, 'One can use either `join_node` or `join_data`'  if options[:join_node] && options[:join_data]
      if options[:join_node]
        raise 'Do not use join_node! It\'s obsolete.'
        join_node = options[:join_node]
      elsif options[:join_data]
        join_node = InnerNode.new(nil, nodes_to_unite, options[:join_data])
      else
        join_node = InnerNode.new(nil, nodes_to_unite, {})
      end

      raise Error, 'One can\'t link non-existent trees'  unless nodes_to_unite.all?{|node| root_nodes.include?(node)}
      Forest.new(@root_nodes - nodes_to_unite + [join_node])
    end

    def +(other)
      case other
      when Forest
        Forest.new(root_nodes.to_a + other.root_nodes.to_a)
      when Tree, Leaf
        Forest.new(root_nodes.to_a + [other])
      end
    end

    def ==(other)
      root_nodes == other.root_nodes
    end

    alias_method :eql?, :==

    def hash
      root_nodes.hash
    end

    def number_of_trees
      root_nodes.size
    end
    def number_of_nodes
      root_nodes.map(&:number_of_nodes).inject(0, &:+)
    end
    def number_of_leafs
      root_nodes.map(&:number_of_leafs).inject(0, &:+)
    end
  end
end
