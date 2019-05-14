require 'clustering_tree'

include Clusterize
describe ClusteringTree do
  #    (0))............(1)
  #    |               |
  # 6  |               |  7
  #    |               |
  #    (2).......\     |
  #               \ ..(3)
  #            8

  let(:max_distance_to_nodes) do
    lambda do |distance_matrix, item_index, node_1_in_cluster_index, node_2_in_cluster_index|
      dist_1 = distance_matrix[item_index][node_1_in_cluster_index]
      dist_2 = distance_matrix[item_index][node_2_in_cluster_index]
      [dist_1, dist_2].max
    end
  end

  let(:item_0) { ClusteringTree::Leaf.new(item: 'leaf 0 - (0,0)', index: 0) }
  let(:item_1) { ClusteringTree::Leaf.new(item: 'leaf 1 - (8,0)', index: 1) }
  let(:item_2) { ClusteringTree::Leaf.new(item: 'leaf 2 - (0,6)', index: 2) }
  let(:item_3) { ClusteringTree::Leaf.new(item: 'leaf 3 - (8,7)', index: 3) }

  let(:distances) do
    [ [0,               8,  6,              Math.sqrt(113)],
      [8,               0,  10,             7],
      [6,               10, 0,              Math.sqrt(65)],
      [Math.sqrt(113),  7,  Math.sqrt(65),  0] ]
  end

  let(:distances_with_one_link) do
    [ [0,               8,  6,              Math.sqrt(113), 6],
      [8,               0,  10,             7,              10],
      [6,               10, 0,              Math.sqrt(65),  6],
      [Math.sqrt(113),  7,  Math.sqrt(65),  0,              Math.sqrt(113)],
      [6,               10, 6,              Math.sqrt(113), 0] ]
  end
  let(:distance_matrix_with_names) { DistanceMatrix.new(distances, ['leaf 0 - (0,0)', 'leaf 1 - (8,0)', 'leaf 2 - (0,6)', 'leaf 3 - (8,7)']) }
  let(:distance_matrix_with_leafs) { DistanceMatrix.new(distances, [item_0, item_1, item_2, item_3]) }
  let(:forest) { Forest.new([item_0, item_1, item_2, item_3]) }

  let(:linkage_node) { Clusterize::Tree.new([item_0, item_2], {item: 'Union (0,2)', index: 4}) }  # Use only as 5-th element.
  let(:distance_matrix_with_one_link) { DistanceMatrix.new(distances_with_one_link, [item_0, item_1, item_2, item_3, linkage_node]) }
  let(:forest_with_one_link) { Forest.new([item_1, item_3, linkage_node]) }

  let(:clustering_tree_without_links) { ClusteringTree.new(distance_matrix_with_leafs, forest) }
  let(:clustering_tree_with_one_link) { ClusteringTree.new(distance_matrix_with_one_link, forest_with_one_link) }

  describe '.new' do
    it 'should raise when distance matrix and forest has different number of items(including inner items)' do
      distance_matrix = DistanceMatrix.new([[0,8],[8,0]], [item_0, item_1])

      flat_forest = Forest.new([item_0, item_1, item_2])
      expect { ClusteringTree.new(distance_matrix, flat_forest) }.to raise_error(InconsistentSizeError)

      same_number_of_roots_forest = Forest.new([item_0, Clusterize::Tree.new([item_1, item_2], 'root')])
      expect { ClusteringTree.new(distance_matrix, same_number_of_roots_forest) }.to raise_error(InconsistentSizeError)
    end

    it 'should not raise when distance matrix and number of nodes in a forest are the same though number of trees is different' do
      union_tree = Clusterize::Tree.new([item_0, item_2], 'union(1,3) - (0,3)')
      distance_matrix = DistanceMatrix.new([[0,8,6,3],[8,0,10, Math.sqrt(73)],[6,10,0,3],[3, Math.sqrt(73), 3, 0]], [item_0, item_1, item_2, union_tree])  # 4 items
      non_flat_forest = Forest.new([item_0, union_tree])  # only 2 trees but 4 nodes (one tree has a root and two leafs)
      expect { ClusteringTree.new(distance_matrix, non_flat_forest) }.not_to raise_error
    end

    context 'with consistent matrix and clustering tree' do
      context 'without links' do
        subject { clustering_tree_without_links }

        its(:distance_matrix) { should eq distance_matrix_with_leafs }
        its(:clustering_forest) { should eq forest }
        its(:items) { should eq [item_0, item_1, item_2, item_3] }
        its(:top_level_items) { should eq Set.new([item_0, item_1, item_2, item_3]) }
      end

      context 'with a link' do
        subject { clustering_tree_with_one_link }

        its(:distance_matrix) { should eq distance_matrix_with_one_link }
        its(:clustering_forest) { should eq forest_with_one_link }
        its(:items) { should eq [item_0, item_1, item_2, item_3, linkage_node] }
        its(:top_level_items) { should eq Set.new([item_1, item_3, linkage_node]) }
      end
    end
  end

  describe '.from_distance_matrix' do
    subject { ClusteringTree.from_distance_matrix(distance_matrix_with_names) }
    its(:distance_matrix) { should eq distance_matrix_with_leafs }
    its(:clustering_forest) { should eq forest }
  end

  describe '#index_of_minimal_element' do
    it 'should return row and column indices for nearest nodes' do
      clustering_tree_without_links.index_of_minimal_element.should == [2,0]
    end
    it 'should return row and column indices for nearest nodes among top level nodes' do
      clustering_tree_with_one_link.index_of_minimal_element.should == [3,1]
    end
  end

  describe '#with_link' do
    let(:tree_w_link) do
      clustering_tree_without_links.with_link(0,2, &max_distance_to_nodes)
    end

    it { tree_w_link.should be_kind_of ClusteringTree }

    it 'has one more item (linkage node added)' do
      tree_w_link.items.should == [item_0, item_1, item_2, item_3, linkage_node]
    end
    it 'has one less top level item (linkage node added, but two items are not on top now)' do
      tree_w_link.top_level_items.should == Set.new([item_1, item_3, linkage_node])
    end
    it 'has distance_matrix expanded' do
      tree_w_link.distance_matrix.should == distance_matrix_with_one_link
    end
  end

  describe '#clusterize' do
    subject { clustering_tree_without_links.clusterize(&max_distance_to_nodes) }
    it { should be_kind_of ClusteringTree }
    it { subject.top_level_items.size.should == 1 }
    it { subject.clustering_forest.number_of_leafs.should == 4 }
  end

  describe '#cluster_tree' do
    subject { clustering_tree_without_links.cluster_tree(&max_distance_to_nodes) }
    it { should be_kind_of ClusteringTree::Tree }
    it { subject.number_of_nodes.should == 2 * 4 - 1 }
    it { subject.number_of_leafs.should == 4 }
  end

  describe '#cluster_tree_with_distances_from_center' do
    subject { clustering_tree_without_links.cluster_tree_with_distances_from_center(&max_distance_to_nodes) }
    it { should be_kind_of ClusteringTree::Tree }
    specify 'each inner node has `distance` data' do
      subject.each_node{|node| node.data.should have_key(:distance)  if !node.leaf? && !node == subject }
    end
    specify 'leaf nodes has no `distance` data' do
      subject.each_leaf{|node| node.data.should_not have_key(:distance) }
    end
    specify 'root node has `distance` data' do
      subject.each_node{|node| node.data.should have_key(:distance) }
    end
  end
end
