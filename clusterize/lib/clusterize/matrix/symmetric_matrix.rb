require_relative '../errors'

module Clusterize
  class SymmetricMatrix
    attr_reader :lower_triangle, :items
    # TODO: identify kind of input: lower_triangle or squared; Possibly create one more factory method for this
    def initialize(lower_triangle, items)
      size = lower_triangle.size
      raise Error, 'Number of items is not the same as matrix size'  unless items.size == size

      raise Error, 'Input is not diagonal lower triangular matrix'  unless self.class.lower_triangle?(lower_triangle)
      raise Error, 'There are duplicate items'  unless items.size == items.uniq.size
      @lower_triangle, @items = lower_triangle, items
    end

    def self.from_square_matrix(matrix, items)
      size = matrix.size
      raise Error, 'Not square matrix' unless matrix.all?{|row| row.size == size}
      raise Error, 'Not symmetric matrix'  unless symmetric?(matrix)
      lower_triangle = matrix.each_with_index.map{|row, row_index| row.first(row_index + 1) }
      self.new(lower_triangle, items)
    end

    def size
      @lower_triangle.size
    end

    def element_at(row, column)
      if row >= column
        lower_triangle[row][column]
      elsif row < column
        lower_triangle[column][row]
      end
    end

    def matrix
      size.times.map do |row_index|
        size.times.map do |column_index|
          element_at(row_index, column_index)
        end
      end
    end

    def [](row, column)
      element_at(row, column)
    end

    def ==(other)
      lower_triangle == other.lower_triangle && items == other.items
    end

    class << self
      def symmetric?(matrix)
        matrix.each_index.all? do |row_index|
          matrix.each_index.all? do |column_index|
            matrix[row_index][column_index] == matrix[column_index][row_index]
          end
        end
      end
      private :symmetric?

      # checks whether array represents lower triangular matrix starting with 1 element vertex (not 0 elements!)
      def lower_triangle?(lower_triangle)
        lower_triangle.each_with_index.all?{|row, row_index|
          row.size == row_index + 1
        }
      end

      # TODO: check whether names are present in a given stream by empty cell in a first row, first column (but make this check able to turn off)
      def from_stream(stream)
        items = stream.readline.chomp.split("\t").drop(1)
        matrix = stream.readlines.map.with_index do |line, index|
          item, *row_elements = line.chomp.split("\t")
          raise Error, 'Items of rows and columns mismatch'  unless item == items[index]
          row_elements.map(&:to_f)
        end
        self.from_square_matrix(matrix, items)
      end

      def from_file(filename)
        File.open(filename) do |filestream|
          from_stream(filestream)
        end
      end
    end
  end
end
