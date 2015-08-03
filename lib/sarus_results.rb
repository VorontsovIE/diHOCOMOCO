# for `besthit suppress` mode
SarusResults = Struct.new(:score, :position, :orientation) do
  def self.from_string(str)
    score, pos, orientation = str.chomp.split("\t")
    self.new(score.to_f, pos.to_i, orientation.to_sym)
  end

  def self.each_in_file(filename)
    return enum_for(:each_in_file, filename)  unless block_given?
    File.open(filename) do |f|
      f.each_line{|line| yield self.from_string(line) }
    end
  end
end
