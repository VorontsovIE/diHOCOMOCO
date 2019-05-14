$:.unshift '../lib'
require 'rspec'
require 'set'
require 'tree'

include Clusterize

describe Tree do
  describe '.new' do
    it 'raises when there are no children' do
      expect { Tree.new([], 'root') }.to raise_error(NoChildrenNodeError)
    end

    it 'raises when there are duplicate children' do
      leaf = Leaf.new('leaf')
      expect { Tree.new([leaf,leaf], 'root') }.to raise_error(DuplicateEntriesError)
    end

    it 'sets data to nil by default' do
      leaf_1 = Leaf.new('leaf_1')
      leaf_2 = Leaf.new('leaf_2')
      Tree.new([leaf_1, leaf_2]).data.should be_nil
    end
  end


  context 'for initialized tree' do
    # leaf(1)    leaf(2)        leaf(3)
    #  \            \-- inner(4) --/
    #   \---- root(5)-----/
    let(:leaf_1) {Leaf.new('data 1')}
    let(:leaf_2) {Leaf.new('data 2')}
    let(:leaf_3) {Leaf.new('data 3')}
    let(:inner_node) { Tree.new([leaf_2, leaf_3], 'inner node') }
    subject { Tree.new([leaf_1, inner_node], 'root node') }

    its(:children) { should include(leaf_1, inner_node) }
    its(:children) { should_not include(leaf_2, leaf_3) }
    its(:data) { should eq 'root node' }
    its(:number_of_nodes) { should eq 5}
    its(:number_of_leafs) { should eq 3}

    its(:subtree_nodes) { should include(leaf_1, leaf_2, leaf_3) }
    its(:subtree_nodes) { should include(inner_node, subject) }
    it{ subject.subtree_nodes.size.should == 5 }

    its(:subtree_leafs) { should include(leaf_1, leaf_2, leaf_3) }
    its(:subtree_leafs) { should_not include(inner_node, subject) }
    it{ subject.subtree_leafs.size.should == 3 }

    it {should_not be_leaf}
    it 'raises on attempt to modify children list' do
      expect { subject.children << Leaf.new('data 4') }.to raise_error
      expect { subject.children = [Leaf.new('data 4')] }.to raise_error
    end
    it 'raises on attempt to modify data' do
      expect { subject.data << ' some other part of data'}.to raise_error
      expect { subject.data = 'another data'}.to raise_error
    end

    describe '#==' do
      it 'should be equal to a tree with the same children and data' do
        subject.should == Tree.new([leaf_1, inner_node], 'root node')
      end
      it 'should be equal to a tree with the same children in different order and the same data' do
        subject.should == Tree.new([inner_node, leaf_1], 'root node')
      end
      it 'should not be equal to a leaf with the same data' do
        subject.should_not == Leaf.new('root node')
      end
      it 'should not be equal to a tree with the different children' do
        subject.should_not == Tree.new([leaf_1, leaf_2], 'root node')
      end
      it 'should not be equal to a tree with different data' do
        subject.should_not == Tree.new([leaf_1, inner_node], 'another root node')
      end
    end

    describe '#eql? and #hash' do
      it 'equal trees are eql' do
        subject.should be_eql Tree.new([leaf_1, inner_node], 'root node')
      end
      it 'not equal forests are not eql' do
        subject.should_not be_eql  Tree.new([leaf_1, leaf_2], 'another root node')
      end
      it 'equal forests have equal hashes' do
        subject.hash.should ==  Tree.new([leaf_1, inner_node], 'root node').hash
      end
      it 'allows work of set operations' do
        Set.new([subject]).should include subject
      end
    end

    describe '#each_node' do
      it 'should enumerate all nodes including inner nodes' do
        enumerated = []
        subject.each_node{|node| enumerated << node}
        enumerated.size.should == 5
        enumerated.should include(leaf_1, leaf_2, leaf_3, inner_node, subject)
      end
    end

    describe '#each_node_with_parent' do
      it 'should enumerate all nodes including inner nodes with their pairs' do
        enumerated = []
        subject.each_node_with_parent{|node, node_parent| enumerated << [node, node_parent] }
        enumerated.size.should == 5
        enumerated.should include([leaf_1, subject], [leaf_2, inner_node], [leaf_3, inner_node], [inner_node, subject], [subject, nil])
      end
      it 'with `parent` option specified substitutes parent for a root' do
        enumerated = []
        subject.each_node_with_parent(parent: 'root node parent'){|node, node_parent| enumerated << [node, node_parent] }
        enumerated.size.should == 5
        enumerated.should include([leaf_1, subject], [leaf_2, inner_node], [leaf_3, inner_node], [inner_node, subject], [subject, 'root node parent'])
      end
    end

    describe '#each_leaf' do
      it 'should enumerate all leaf-nodes but not inner nodes' do
        enumerated = []
        subject.each_leaf{|node| enumerated << node}
        enumerated.size.should == 3
        enumerated.should include(leaf_1, leaf_2, leaf_3)
        enumerated.should_not include(inner_node, subject)
      end
    end

    describe '#map_tree' do
      it 'should create new tree with the same structure but different data (obtained from block call on nodes)' do
        inner_node_with_suffix = Tree.new([Leaf.new('data 2 with suffix'), Leaf.new('data 3 with suffix')], 'inner node with suffix')
        subject.map_tree{|node| node.data + ' with suffix' }.should == Tree.new([Leaf.new('data 1 with suffix'), inner_node_with_suffix], 'root node with suffix')
      end
    end

    # describe '#subtree_clusters' do
    #   it { subject.subtree_clusters {|node| node.data == 'inner_node' }.should == Forest.new([leaf_1, leaf_2, leaf_3]) }
    #   it { subject.subtree_clusters {|node| node.data == 'root node' }.should == Forest.new([leaf_1, inner_node]) }
    #   it { subject.subtree_clusters {|node| node.data == 'leaf_1' }.should == Forest.new([leaf_1, leaf_2, leaf_3]) } ## not yet working
    #   it { subject.subtree_clusters {|node| false }.should == Forest.new([subject]) }
    # end
  end
end

describe Leaf do
  context 'for initialized leaf' do
    # leaf(1)    leaf(2)        leaf(3)
    #  \            \-- inner(4) --/
    #   \---- root(5)-----/
    subject { Leaf.new('leaf data') }

    its(:children) { should eq [] }
    its(:data) { should eq 'leaf data' }
    its(:number_of_nodes) { should eq 1}
    its(:number_of_leafs) { should eq 1}

    its(:subtree_nodes) { should eq [subject] }
    its(:subtree_leafs) { should eq [subject] }

    it {should be_leaf}

    it 'raises on attempt to modify data' do
      expect { subject.data << ' some other part of data'}.to raise_error
      expect { subject.data = 'another data'}.to raise_error
    end

    describe '#==' do
      it 'should be equal to a leaf with the same data' do
        subject.should == Leaf.new('leaf data')
      end
      it 'should not be equal to a leaf with different data' do
        subject.should_not == Leaf.new('another node')
      end
      it 'should not be equal to tree with the same data' do
        subject.should_not == Tree.new([Leaf.new('some child')], 'leaf data')
      end
    end

    describe '#eql? and #hash' do
      it 'equal leafs are eql' do
        subject.should be_eql Leaf.new('leaf data')
      end
      it 'not equal forests are not eql' do
        subject.should_not be_eql Leaf.new('another leaf data')
      end
      it 'equal forests have equal hashes' do
        subject.hash.should == Leaf.new('leaf data').hash
      end
      it 'allows work of set operations' do
        Set.new([subject]).should include subject
      end
    end

    describe '#each_node' do
      it 'should yield itself' do
        expect{|b| subject.each_leaf(&b) }.to yield_with_args(subject)
      end
    end
    describe '#each_leaf' do
      it 'should yield itself' do
        expect{|b| subject.each_leaf(&b) }.to yield_with_args(subject)
      end
    end

    describe '#map_tree' do
      it 'should create new leaf with data from called block' do
        subject.map_tree{|node| node.data + ' with suffix' }.should == Leaf.new('leaf data with suffix')
      end
    end

    # describe '#subtree_clusters' do
    # end
  end
end

describe Forest do
  #               |                  |
  #             root(4)            leaf(3)
  #  leaf(1)  --/    \--leaf(2)
  let(:leaf_1) { Leaf.new('leaf 1') }
  let(:leaf_2) { Leaf.new('leaf 2') }
  let(:leaf_3) { Leaf.new('leaf 3') }
  let(:leaf_4) { Leaf.new('leaf 4') }
  let(:tree_1) { Tree.new([leaf_1, leaf_2], 'tree 1') }
  let(:tree_2) { leaf_3 }
  let(:tree_3) { leaf_4 }

  subject { Forest.new([tree_1, tree_2]) }
  its(:trees) {should include(tree_1, tree_2) }
  its(:number_of_nodes) {should eq 4}
  its(:number_of_trees) {should eq 2}
  its(:number_of_leafs) {should eq 3}

  describe '.new ' do
    it 'should raise when there are duplicate trees' do
      expect { Forest.new([tree_1, tree_2, tree_1]) }.to raise_error(DuplicateEntriesError)
    end
  end

  describe '#==' do
    it 'should be equal to a forest with the same trees' do
      subject.should == Forest.new([tree_1, tree_2])
    end
    it 'should be equal to a forest with the same trees in different order' do
      subject.should == Forest.new([tree_2, tree_1])
    end
    it 'should not be equal to a forest with the same trees' do
      subject.should_not == Forest.new([tree_3, tree_2])
    end
  end

  describe '#eql? and #hash' do
    it 'equal forests are eql' do
      Forest.new([tree_1, tree_2]).should be_eql Forest.new([tree_1, tree_2])
    end
    it 'not equal forests are not eql' do
      Forest.new([tree_1, tree_2]).should_not be_eql Forest.new([tree_1, tree_3])
    end
    it 'equal forests have equal hashes' do
      Forest.new([tree_1, tree_2]).hash.should == Forest.new([tree_2, tree_1]).hash
    end
    it 'allows work of set operations' do
      Set.new([subject]).should include subject
    end
  end

  describe '#+' do
    context 'with forest argument' do
      it 'creates a tree consisting of all trees of both forests' do
        sum = subject + Forest.new([tree_3])
        sum.should == Forest.new([tree_1, tree_2, tree_3])
      end
      it 'raises when summand forests have the same trees' do
        expect{ subject + Forest.new([tree_2, tree_3]) }.to raise_error(DuplicateEntriesError)
      end
      it 'raises when summand tree is already in the forest' do
        expect{ subject + tree_1 }.to raise_error(DuplicateEntriesError)
        expect{ Forest.new([tree_1, leaf_3]) + leaf_3 }.to raise_error(DuplicateEntriesError)
      end
    end
    context 'with Tree argument' do
      it 'returns a forest with the tree added' do
        sum = subject + tree_3
        sum.should == Forest.new([tree_1, tree_2, tree_3])
      end
    end
    context 'with leaf-Tree argument' do
      it 'returns a forest with the tree added' do
        sum = subject + leaf_4
        sum.should == Forest.new([tree_1, tree_2, leaf_4])
      end
    end
  end


  describe '#with_link' do
    it 'unites several trees in a single tree' do
      new_forest = subject.with_link([tree_1, tree_2])
      new_forest.trees.should_not include(tree_1, tree_2)
      new_forest.trees.should include Tree.new([tree_1,tree_2])
    end

    it 'without `join_data` or `join_node` options joined tree has `nil` data' do
      subject.with_link([tree_1, tree_2]).trees.should include Tree.new([tree_1,tree_2], nil)
    end
    it 'with `join_data` option, join-tree has specified data' do
      subject.with_link([tree_1, tree_2], join_data: 'union data').trees.should include Tree.new([tree_1,tree_2], 'union data')
    end
    it 'with `join_node` option, join-tree is specified by this option' do
      join_node = Tree.new([tree_1, tree_2], 'union data')
      subject.with_link([tree_1, tree_2], join_node: join_node).trees.should include join_node
    end
    it 'with both `join_node` and `join_data` options it raises' do
      join_node = Tree.new([tree_1, tree_2], 'union data')
      expect{ subject.with_link([tree_1, tree_2], join_node: join_node, join_data: 'union data') }.to raise_error(IllegalOptionsError)
    end

    it 'unites only given trees into a single node' do
      new_forest = Forest.new([tree_1, tree_2, tree_3]).with_link([tree_1, tree_2])
      new_forest.trees.should_not include(tree_1, tree_2)
      new_forest.trees.should include Tree.new([tree_1,tree_2])
      new_forest.trees.should include tree_3
    end
    it 'raises if nodes to unite don\'t exist in a forest' do
      expect { subject.with_link([tree_1, tree_3]) }.to raise_error NoSuchNodeError
    end
  end
end
