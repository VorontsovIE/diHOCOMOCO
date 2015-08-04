# threshold, pvalue pairs
# num_words is printed only to terminal
FindPvalueResults = Struct.new(:threshold, :pvalue) do
  def self.from_string(str)
    threshold, pvalue = str.chomp.split("\t")
    self.new(threshold.to_f, pvalue.to_f)
  end

  def self.each_in_file(filename)
    return enum_for(:each_in_file, filename)  unless block_given?
    File.open(filename) do |f|
      f.each_line.reject{|line| line.start_with?('#') }.each{|line| yield self.from_string(line) }
    end
  end
end
