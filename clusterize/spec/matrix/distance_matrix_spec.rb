require 'stringio'
require 'clusterize/matrix/distance_matrix'

include Clusterize

describe DistanceMatrix do
  describe '.from_square_matrix' do
    it 'should raise on non-symmetric matrix' do
      expect{ DistanceMatrix.from_square_matrix([[0,2],[3,0]], ['first item', 'second item']) }.to raise_error(Clusterize::Error)
    end
    it 'should raise if matrix diagonal is not zero' do
      expect{ DistanceMatrix.from_square_matrix([[0,2],[2,1]], ['first item', 'second item']) }.to raise_error(Clusterize::Error)
    end
    it 'should raise if there are negative matrix elements' do
      expect{ DistanceMatrix.from_square_matrix([[0,-2],[-2,0]], ['first item', 'second item']) }.to raise_error(Clusterize::Error)
    end
  end
  describe '.from_stream' do
    it 'should return DistanceMatrix' do
      stream = StringIO.new("\tfirst\tsecond\n" +
                            "first\t0\t2\n" +
                            "second\t2\t0\n")
      expect(DistanceMatrix.from_stream(stream)).to be_kind_of DistanceMatrix
    end
  end

  context 'valid matrix' do
    let(:matrix) { [[0,2,4],[2,0,1],[4,1,0]] }
    let(:items) { ['first item', 'second item', 'third item'] }
    subject { DistanceMatrix.from_square_matrix(matrix, items) }

    specify { expect(subject.to_distance_matrix).to eq subject }

    it 'should be equal to distance matrix with the same matrix and items' do
      expect(subject).to eq DistanceMatrix.from_square_matrix([[0,2,4],[2,0,1],[4,1,0]], ['first item', 'second item', 'third item'])
    end
    it 'should not be equal to common-type matrix' do
      expect(subject).not_to eq SymmetricMatrix.from_square_matrix([[0,2,4],[2,0,1],[4,1,0]], ['first item', 'second item', 'third item'])
    end

    describe '#index_of_minimal_element' do
      it 'should return array with two indices' do
        expect(subject.index_of_minimal_element).to be_kind_of Array
        expect(subject.index_of_minimal_element.size).to eq 2
      end
      it 'should not be on diagonal' do
        row_index, column_index = subject.index_of_minimal_element
        expect(row_index).not_to eq column_index
      end
      it 'should be behind diagonal' do
        row_index, column_index = subject.index_of_minimal_element
        expect(row_index).to be > column_index
      end
      it 'should give index of lowest element' do
        row_index, column_index = subject.index_of_minimal_element
        expect(row_index).to eq 2
        expect(column_index).to eq 1
      end
      it 'when indices subset is given, it choose minimum amongst items in a subset' do
        row_index, column_index = subject.index_of_minimal_element([0,2])
        expect(row_index).to eq 2
        expect(column_index).to eq 0
      end
    end

    describe '#with_node' do
      it 'fails if size of row is not equal to matrix size' do
        expect{ subject.with_node([1,2], 'fourth item') }.to raise_error(Clusterize::Error)
        expect{ subject.with_node([1,2,3,4], 'fourth item') }.to raise_error(Clusterize::Error)
      end
      it 'inserts new item as last item and sets distances to already present items' do
        expect(subject.with_node([10,20,30], 'fourth item')).to eq DistanceMatrix.from_square_matrix([[0,2,4, 10],[2,0,1, 20],[4,1,0, 30],[10, 20, 30, 0]], ['first item', 'second item', 'third item', 'fourth item'] )
      end
    end
  end

  describe '#index_of_minimal_element' do
    let(:distance_matrix) { DistanceMatrix.from_square_matrix([[0]], ['the only item']) }
    context 'for single-element matrix' do
      it 'raises an exception' do
        expect { distance_matrix.index_of_minimal_element }.to raise_error(Clusterize::Error)
      end
    end
  end
end
