class IUPAC < String
  CODE = {"A" => "A", "C" => "C", "G" => "G", "T" => "T", 
          "AG" => "R", "CT" => "Y", "GT" => "K", "AC" => "M", 
          "CG" => "S", "AT" => "W", "CGT" => "B", "AGT" => "D", "ACT" => "H", "ACG" => "V", "ACGT" => "N"}
  REVCODE = CODE.invert
  
  def dup
    IUPAC.new(self)
  end
  
  def initialize(words)
    if words.is_a?(Array)
      iupac = (0...words[0].size).collect { |i|
        (0...words.size).collect { |j| words[j][i,1] }.uniq.sort.inject("") { |cola, letter| cola += letter }
      }.inject("") { |iup, cola|
        checkerr("bad letter set #{cola}") { !CODE.has_key?(cola) }
        iup += CODE[cola]
      }
      super(iupac)
    elsif words.is_a?(IUPAC)
      super(words)
    elsif words.is_a?(String)
      checkerr("word #{words} has strange characters") { words.tr('ACGTURYKMSWBDHVN', '').size > 0 }
      super(words)
    end
  end
  
  def ==(iupac)
    return false if self.size != iupac.size
    (0...self.size).inject(true) { |result, i| result &= IUPACOM[self[i,1]][iupac[i,1]] }
  end
  
  def merge(iupac)
    return nil if self.size != iupac.size
    res = (0...self.size).inject("") { |res, i|
      merges = REVCODE[self[i,1]].split(//).concat(REVCODE[iupac[i,1]].split(//)).uniq.sort.inject("") { |s, c| s += c}
      res << CODE[merges]
    }
    return IUPAC.new(res)
  end
  
  def include?(iupac)
    return false if self.size < iupac.size || !iupac.is_a?(IUPAC)
    (0..self.size-iupac.size).each { |i|
      return i if IUPAC.new(self[i,iupac.size]) == iupac
    }
    return false
  end
  
  def compl
    return self.tr("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")
  end
  
  def compl!
    self.tr!("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")
    return self
  end
  
  alias reverse_string reverse
  def reverse
    return IUPAC.new(reverse_string)
  end
  
  alias comp! compl!
  alias complement! compl!
  alias comp compl
  alias complement compl
  
private
  IUPACOM = { "A" => {"A" => :llib, "R" => :llib, "M" => :llib, "W" => :llib, "D" => :llib, "H" => :llib, "V" => :llib, "N" => :llib},
                "C" => {"C" => :llib, "Y" => :llib, "M" => :llib, "S" => :llib, "B" => :llib, "H" => :llib, "V" => :llib, "N" => :llib},
                "G" => {"G" => :llib, "R" => :llib, "K" => :llib, "S" => :llib, "B" => :llib, "D" => :llib, "V" => :llib, "N" => :llib},
                "T" => {"T" => :llib, "Y" => :llib, "K" => :llib, "W" => :llib, "B" => :llib, "D" => :llib, "H" => :llib, "N" => :llib}
  }
  IUPACOM["R"] = IUPACOM["G"].merge(IUPACOM["A"])
  IUPACOM["Y"] = IUPACOM["T"].merge(IUPACOM["C"])
  IUPACOM["K"] = IUPACOM["G"].merge(IUPACOM["T"])
  IUPACOM["M"] = IUPACOM["A"].merge(IUPACOM["C"])
  IUPACOM["S"] = IUPACOM["G"].merge(IUPACOM["C"])
  IUPACOM["W"] = IUPACOM["A"].merge(IUPACOM["T"])
  IUPACOM["B"] = IUPACOM["G"].merge(IUPACOM["T"].merge(IUPACOM["C"]))
  IUPACOM["D"] = IUPACOM["G"].merge(IUPACOM["A"].merge(IUPACOM["T"]))
  IUPACOM["H"] = IUPACOM["A"].merge(IUPACOM["C"].merge(IUPACOM["T"]))
  IUPACOM["V"] = IUPACOM["G"].merge(IUPACOM["C"].merge(IUPACOM["A"]))
  IUPACOM["N"] = IUPACOM["A"].merge(IUPACOM["C"].merge(IUPACOM["G"].merge(IUPACOM["T"])))
 
#  IUPACMERGE = CODE.merge({
#    "AA" => "A", "CC" => "C", "GG" => "G", "TT" => "T",
#    
#  })

end
