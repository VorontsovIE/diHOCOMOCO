require 'stringio'
require 'clusterize/matrix/distance_matrix'

include Clusterize
describe SymmetricMatrix do
  describe '.new' do
    let(:items) {  }
    specify 'raises on non-triangular matrices' do  # Probably will be changed to automatically guess input
      expect{SymmetricMatrix.new([[1,2,4],
                                  [2,3,5],
                                  [4,5,6]], ['first item', 'second item', 'third item']) }.to raise_error(Clusterize::Error)
    end
    specify 'raises on matrices with zero elements vertex' do
      expect{SymmetricMatrix.new([[],
                                  [1],
                                  [2,3]], ['first item', 'second item', 'third item']) }.to raise_error(Clusterize::Error)
    end
    specify 'raises on non-triangular matrices' do
      expect{SymmetricMatrix.new([[1],
                                  [2,3],
                                  [4,5]], ['first item', 'second item', 'third item']) }.to raise_error(Clusterize::Error)
    end
    specify 'raises on upper-triangular matrices' do  # Probably will be changed to automatically guess input
      expect{SymmetricMatrix.new([[1,2,4],
                                  [  3,5],
                                  [    6]], ['first item', 'second item', 'third item']) }.to raise_error(Clusterize::Error)
    end
    specify 'raises when number of items not equal to matrix size' do
      expect{SymmetricMatrix.new([[1],
                                  [2,3]], ['first item', 'second item', 'third item']) }.to raise_error(Clusterize::Error)
    end
    specify 'raises when items are duplicated' do
      expect{SymmetricMatrix.new([[1],
                                  [2,3],
                                  [4,5,6]], ['first item', 'second item', 'second item']) }.to raise_error(Clusterize::Error)
    end

    context 'correct matrix' do
      let(:lower_triangle) { [[1],
                              [2,3],
                              [4,5,6]] }
      let(:items) { ['first item', 'second item', 'third item'] }
      let(:matrix) { SymmetricMatrix.new(lower_triangle, items)}
      specify 'matrix can be initialized' do
        expect{ matrix }.not_to raise_error
      end
      specify{ expect(matrix.lower_triangle).to eq lower_triangle }
      specify{ expect(matrix.items).to eq items }
      specify{
        expect(matrix.element_at(0,0)).to eq 1
        expect(matrix.element_at(0,1)).to eq 2
        expect(matrix.element_at(1,0)).to eq 2
      }

      specify {
        expect(matrix.matrix).to eq [ [1,2,4],
                                      [2,3,5],
                                      [4,5,6] ]

      }
    end

  end
  describe '.from_square_matrix' do
    it 'should raise on rectangular(non-square matrices)' do
      expect{ SymmetricMatrix.from_square_matrix([[1,2],[3,4],[5,6]], ['first item', 'second item']) }.to raise_error(Clusterize::Error)
      expect{ SymmetricMatrix.from_square_matrix([[1,2,3],[4,5,6]], ['first item', 'second item']) }.to raise_error(Clusterize::Error)
    end
    it 'should raise on non-matrix array (with different row sizes)' do
      expect{ SymmetricMatrix.from_square_matrix([[1,2],[3,4,5]], ['first item', 'second item']) }.to raise_error(Clusterize::Error)
    end
    it 'should raise if items list has size different with size of matrix' do
      expect{ SymmetricMatrix.from_square_matrix([[1,2],[4,5]], ['first item', 'second item', 'third item']) }.to raise_error(Clusterize::Error)
    end
    it 'should raise if there are duplicate items' do
      expect{ SymmetricMatrix.from_square_matrix([[0,2,2],[2,0,0], [2,0,0]], ['first item', 'second item', 'second item']) }.to raise_error(Clusterize::Error)
    end
    it 'should raise if matrix is not symmetric' do
      expect{ SymmetricMatrix.from_square_matrix([[0,2,2],[2,0,0], [2,10,0]], ['first item', 'second item', 'third item']) }.to raise_error(Clusterize::Error)
    end
  end

  describe '.from_stream' do
    it 'should parse a matrix' do
      stream = StringIO.new("\tfirst\tsecond\n" +
                            "first\t1\t2\n" +
                            "second\t2\t1\n")
      matrix = SymmetricMatrix.from_stream(stream)
      expect(matrix.matrix).to eq [[1,2],[2,1]]
      expect(matrix.lower_triangle).to eq [[1],[2,1]]
      expect(matrix.items).to eq ['first', 'second']
    end
    it 'should parse matrix with spaces in node names' do
      stream = StringIO.new("\tfirst item\tsecond item\n" +
                            "first item\t1\t2\n" +
                            "second item\t2\t1\n")
      matrix = SymmetricMatrix.from_stream(stream)
      expect(matrix.matrix).to eq [[1,2],[2,1]]
      expect(matrix.lower_triangle).to eq [[1],[2,1]]
      expect(matrix.items).to eq ['first item', 'second item']
    end
    it 'should raise if row and column nodes are different' do
      stream = StringIO.new("\tfirst\tsecond\n" +
                            "third\t1\t2\n" +
                            "fourth\t2\t1\n")
      expect{ SymmetricMatrix.from_stream(stream) }.to raise_error(Clusterize::Error)
    end
    it 'should raise if row and column nodes are in different order', :to_be_changed => true do
      stream = StringIO.new("\tfirst\tsecond\n" +
                            "second\t1\t2\n" +
                            "first\t2\t1\n")
      expect{ SymmetricMatrix.from_stream(stream) }.to raise_error(Clusterize::Error)
    end

  end

  context 'valid matrix' do
    let(:matrix) { [[1,2,3],
                    [2,1,4],
                    [3,4,1.5]] }
    let(:different_matrix) { [[1,2,3],
                              [2,100,4],
                              [3,4,1.5]] }
    let(:items) { ['first item', 'second item', 'third item'] }
    let(:different_items) { ['first item', 'second item', 'another item'] }
    subject { SymmetricMatrix.from_square_matrix(matrix, items) }

    specify{ expect(subject.lower_triangle).to eq [ [1],
                                                    [2, 1],
                                                    [3, 4, 1.5] ]
    }
    specify{ expect(subject.matrix).to eq matrix }
    specify{ expect(subject.items).to eq items }
    specify{ expect(subject.size).to eq 3 }
    specify '#element_at(row index,column index) returns element (indices 0-based)' do
      expect(subject.element_at(0,0)).to eq 1
      expect(subject.element_at(0,1)).to eq 2
      expect(subject.element_at(0,2)).to eq 3

      expect(subject.element_at(1,0)).to eq 2
      expect(subject.element_at(1,1)).to eq 1
      expect(subject.element_at(1,2)).to eq 4

      expect(subject.element_at(2,0)).to eq 3
      expect(subject.element_at(2,1)).to eq 4
      expect(subject.element_at(2,2)).to eq 1.5
    end
    specify '#[row index,column index] returns element (indices 0-based)' do
      expect(subject[0,0]).to eq 1
      expect(subject[0,1]).to eq 2
      expect(subject[0,2]).to eq 3

      expect(subject[1,0]).to eq 2
      expect(subject[1,1]).to eq 1
      expect(subject[1,2]).to eq 4

      expect(subject[2,0]).to eq 3
      expect(subject[2,1]).to eq 4
      expect(subject[2,2]).to eq 1.5
    end

    it 'should be equal to matrix with the same matrix and items' do
      expect(subject).to eq SymmetricMatrix.from_square_matrix(matrix, items)
    end
    it 'should not be equal to matrix with the different matrix' do
      expect(subject).not_to eq SymmetricMatrix.from_square_matrix(different_matrix, items)
    end
    it 'should not be equal to matrix with the different items' do
      expect(subject).not_to eq SymmetricMatrix.from_square_matrix(matrix, different_items)
    end
  end
end
