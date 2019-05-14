module Tree
  module TreeNode
    # Tag module
    class << self
      protected :new
      def create(options = {})
        if options.has_key?(:children) && !options[:children].empty?
          InnerNode.new(options[:children], options.fetch(:data, {}))
        else
          LeafNode.new(options.fetch(:data, {}))
        end
      end
    end

  end
end
