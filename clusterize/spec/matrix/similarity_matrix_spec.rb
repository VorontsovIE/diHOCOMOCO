require 'stringio'
require 'clusterize/matrix/similarity_matrix'

include Clusterize
describe SimilarityMatrix do
  describe '.from_square_matrix' do
    it 'should raise on non-symmetric matrix' do
      expect{ SimilarityMatrix.from_square_matrix([[1.0, 0.2],[0.3, 1.0]], ['first item', 'second item']) }.to raise_error(Clusterize::Error)
    end
    it 'should raise if matrix diagonal is not unity' do
      expect{ SimilarityMatrix.from_square_matrix([[0.8, 0.2],[0.2, 1.0]], ['first item', 'second item']) }.to raise_error(Clusterize::Error)
    end
    it 'should raise if there are negative matrix elements' do
      expect{ SimilarityMatrix.from_square_matrix([[1.0, -0.2],[-0.2, 1.0]], ['first item', 'second item']) }.to raise_error(Clusterize::Error)
    end
    it 'should raise if there are matrix elements greater than 1' do
      expect{ SimilarityMatrix.from_square_matrix([[1.0, 1.5],[1.5, 1.0]], ['first item', 'second item']) }.to raise_error(Clusterize::Error)
    end
  end

  describe '.from_stream' do
    it 'should return SimilarityMatrix' do
      stream = StringIO.new("\tfirst\tsecond\n" +
                            "first\t1\t0.6\n" +
                            "second\t0.6\t1\n")
      expect(SimilarityMatrix.from_stream(stream)).to be_kind_of SimilarityMatrix
    end
  end

  context 'valid matrix' do
    let(:matrix) { [[1.0, 0.2],[0.2, 1.0]] }
    let(:items) { ['first item', 'second item'] }
    subject { SimilarityMatrix.from_square_matrix(matrix, items) }

    specify{ expect(subject.to_distance_matrix).to eq DistanceMatrix.from_square_matrix([[0, 0.8],[0.8,0]], items) }

    it 'should be equal to similarity matrix with the same matrix and items' do
      expect(subject).to eq SimilarityMatrix.from_square_matrix([[1.0, 0.2],[0.2, 1.0]] , ['first item', 'second item'])
    end
    it 'should not be equal to common-type matrix' do
      expect(subject).not_to eq SymmetricMatrix.from_square_matrix([[1.0, 0.2],[0.2, 1.0]] , ['first item', 'second item'])
    end
  end
end
