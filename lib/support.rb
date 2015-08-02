Nucleotides = %w[A C G T]
Dinucleotides = Nucleotides.product(Nucleotides).map(&:join)

module InitializeFromHash
  def normalized_hash(hsh)
    result = hsh.map{|k,v|
      [k.to_s.upcase, v]
    }.to_h
    result.default = 0.0
    result
  end

  private :normalized_hash
end
