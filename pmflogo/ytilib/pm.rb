#!/usr/bin/ruby
module Ytilib

require 'ytilib.rb'

#report "PM#save supports simplified .pat format (.pwm as .pat without header)", "pm"
#report "PM#score works only with CAPITAL letters in words", "pm"
#report "TODO: implement basic small-BiSMark validation", "pm"
#report "TODO: PM.load (and new) incorrectly assigns words_count for PWMs", "pm"
#report "INFO: PPM.iupacomp! works SIMILARLY to PM.iupacomp! (so ppm.score finally works correctly)"
#report "TODO: need to include background parameter for the from_IUPAC functions"

class PM

  attr_reader :matrix, :size
  attr_accessor :words_count
  
  alias length size
  
  def score_mean(bckgr = Randoom::DEF_PROBS)
    (0...@size).inject(0.0) { |mean, i| mean += ['A','C','G','T'].inject(0.0) { |sum,l| sum += @matrix[l][i] * bckgr[l] } }
  end
  
  def score_variance(bckgr = Randoom::DEF_PROBS)
    (0...@size).inject(0.0) { |m2, i| 
      deltai = ['A','C','G','T'].inject(0.0) { |sum,l| sum += @matrix[l][i]**2 * bckgr[l] } - ['A','C','G','T'].inject(0.0) { |sum,l| sum += matrix[l][i] * bckgr[l] }**2
      m2 += deltai
    }
  end
  
  def p_value(threshold, mean = nil, variance = nil)
    mean = mean ? mean : score_mean
    variance = variance ? variance : score_variance
    n_ = (threshold - mean) / Math.sqrt(variance)
    p_value = (1 - Math.erf2(n_/Math.sqrt(2))) / 2.0
  end
  
  def best_word
    return (0...size).inject("") { |word, i|
      max = ['A', 'C', 'G', 'T'].collect { |l| @matrix[l][i] }.max
      maxlets = ['A', 'C', 'G', 'T'].select { |l| @matrix[l][i] == max }
      word << (maxlets.size == 1 ? maxlets.first : "N")
    }
  end
  
  def strict_consensus
    return IUPAC.new((0...size).inject("") { |word, i|
      max = ['A', 'C', 'G', 'T'].collect { |l| @matrix[l][i] }.max
      maxlets = ['A', 'C', 'G', 'T'].inject("") { |lets, l| lets += @matrix[l][i] == max ? l : ""}
      word += IUPAC::CODE[maxlets]
    })
  end
  
  def consensus_string(beautiful = false)
    checkerr("words count is undefined") { !@words_count }
    i2o4, thc, tlc = icd2of4, icdThc, icdTlc
    icd = infocod
    
    return String.new((0...size).inject("") { |word, i|

      scores = ['A', 'C', 'G', 'T'].collect { |l| @matrix[l][i] }.uniq.sort.reverse
      
      if icd[i] > i2o4
        scores = [scores.first]
      elsif icd[i] > thc
        scores = scores[0..1]
      elsif icd[i] > tlc
        scores = scores[0..2]
      end
      
      lets = ['A', 'C', 'G', 'T'].inject("") { |lets, l| lets += scores.include?(@matrix[l][i]) ? l : ""}
      
      reslet = IUPAC::CODE[lets]
      reslet = reslet.downcase if beautiful && lets.size > 2
      
      word += reslet
    })
  end
  
  def consensus
    checkerr("words count is undefined") { !@words_count }
    i2o4, thc, tlc = icd2of4, icdThc, icdTlc
    icd = infocod
    
    return IUPAC.new((0...size).inject("") { |word, i|

      scores = ['A', 'C', 'G', 'T'].collect { |l| @matrix[l][i] }.uniq.sort.reverse
      
      if icd[i] > i2o4
        scores = [scores.first]
      elsif icd[i] > thc
        scores = scores[0..1]
      elsif icd[i] > tlc
        scores = scores[0..2]
      end
      
      lets = ['A', 'C', 'G', 'T'].inject("") { |lets, l| lets += scores.include?(@matrix[l][i]) ? l : ""}
      
      word += IUPAC::CODE[lets]
    })
  end
  
  def find_hit(s, score_g, use2strands = true)
    (0..(s.size - @size)).each { |i|
      seq, seq_rc = s[i, @size], s[i, @size].revcomp!
      score_p, score_rc = score(seq), score(seq_rc)
      r = use2strands ? [score_p,score_rc].max : score_p
      return i if r >= score_g
    }
    return nil
  end
  
  def find_hits(s, score_g, use2strands = true)
    (0..(s.size - @size)).select { |i|
      seq, seq_rc = s[i, @size], s[i, @size].revcomp!
      score_p, score_rc = score(seq), score(seq_rc)
      r = use2strands ? [score_p,score_rc].max : score_p
      r >= score_g ? i : nil
    }.compact
  end
  
  def collect_hits(s, score_g, use2strands = true)
    result = []
    (0..(s.size - @size)).each { |i|
      seq, seq_rc = s[i, @size], s[i, @size].revcomp!
      score_p, score_rc = score(seq.upcase), score(seq_rc.upcase)
      result << [score_p, seq, false, i] if score_p >= score_g
      result << [score_rc, seq_rc, true, i] if score_rc >= score_g
    }
    result
  end
  
  def best_hit(s, use2strands = true)
    
    checkerr("too short sequence") { s.size < @size }
    return (0..(s.size - @size)).inject(-Float::MAX) { |r, i|
      seq, seq_rc = s[i, @size], s[i, @size].revcomp!
      score_p, score_rc = score(seq), score(seq_rc)
      r = use2strands ? [r,score_p,score_rc].max : [r,score_p].max
    }
  end
  
  def eql?(pm)
    return ['A','C','G','T'].inject(true) { |equal, letter|
      equal = equal && @matrix[letter].eql?(pm.matrix[letter])
    }
  end
  
  def flexeql?(pm)
    checkerr("for what?") { true }
    return ['A','C','G','T'].inject(true) { |equal, letter|
      # report "letter=#{letter}"
      equal = equal && (0...@size).inject(true) { |deepequal, position| 
        # report "position=#{position}, delta=#{@matrix[letter][position] - pm.matrix[letter][position]}"
        deepequal = deepequal && (@matrix[letter][position] - pm.matrix[letter][position]).abs < 10**-11 
      }
    }
  end
  
  def initialize(size, matrix = nil, words_count = nil)
    checkerr("matrix['A'].size != size, #{matrix['A'].size} != #{size}") { matrix != nil && size != matrix['A'].size }
    @size = size
    @matrix = matrix == nil ? PM.new_matrix(size) : matrix
    if !words_count || words_count <= 0
      words_count = col_sum(0)
      @words_count = words_count.round >= 2 ? words_count.round : nil
    else
      @words_count = words_count
    end
  end
  
  def col_sum(index = 0, letset = ['A','C','G','T'])
    return letset.inject(0) { |sum, l| sum += @matrix[l][index] }
  end
  
  def PM.col_sum(matrix, index = 0)
    return matrix['A'][index] + matrix['C'][index] + matrix['G'][index] + matrix['T'][index]
  end
  	
  def to_pwm!(words_count = nil, probs = Randoom::DEF_PROBS, pseudocount = 1)
    @words_count = words_count if words_count && words_count > 0
    
    @matrix.each_key do |letter|
      (0...@size).each { |pos|
        
        #p "pcm"
        #p @matrix[letter][pos]
        #p @matrix[letter][pos] + (probs[letter] * pseudocount)
        #p ( (@words_count + pseudocount) * probs[letter])
        #exit
        
        @matrix[letter][pos] = Math::log( (@matrix[letter][pos] + (probs[letter] * pseudocount)) / ( (@words_count + pseudocount) * probs[letter]) )
        
      }
    end
    
    return self
  end
  
  def get_pwm(words_count = nil, probs = Randoom::DEF_PROBS, pseudocount = 1)
    return self.dup.to_pwm!(words_count, probs, pseudocount)
  end
  alias to_pwm get_pwm
  
  def get_ppm(words_count = nil)
    words_count = @words_count unless words_count
    checkerr("undefined words count") { !words_count || words_count <= 0 }
    ppm = @matrix['N'] ?  PM.new_matrix_iupac(@size) : PM.new_matrix(@size)
    @matrix.each_key { |letter|
      (0...@size).each { |i|
        ppm[letter][i] = @matrix[letter][i].to_f / words_count
      }
    }
    return PPM.new(@size, ppm, words_count)
  end
  alias to_ppm get_ppm
  
  def get_safe_ppm(words_count = nil) # for matrices w. unbalanced word counts
    words_count = @words_count unless words_count
    checkerr("undefined words count") { !words_count || words_count <= 0 }
    ppm = @matrix['N'] ?  PM.new_matrix_iupac(@size) : PM.new_matrix(@size)
    @matrix.each_key { |letter|
      (0...@size).each { |i|
        sumi = ['A','C','G','T'].collect { |l| @matrix[l][i] }.sum
        ppm[letter][i] = @matrix[letter][i].to_f / sumi
      }
    }
    return PPM.new(@size, ppm, words_count)
  end
  
  def score(word)
    checkerr("word size != pwm.size") { @size != word.size }
    checkerr("word #{word} has strange characters") { 
      @matrix.keys.include?('N') ? word.tr('ACGTRYKMSWBDHVN', '').size > 0 : word.tr('ACGT', '').size > 0
    }
    return (0...@size).inject(0) { |sum, i| 
      sum += @matrix[word[i,1]][i]
    }
  end
	
  def best_score
    return (0...size).inject(0) { |sum, i|
      sum += ['A', 'C', 'G', 'T'].collect { |l| @matrix[l][i] }.max
    }
  end
	
  def worst_score
    return (0...size).inject(0) { |sum, i|
      sum += ['A', 'C', 'G', 'T'].collect { |l| @matrix[l][i] }.min
    }
  end
  
  def dup
    new_matrix = {}
    @matrix.each_key { |letter| new_matrix[letter] = @matrix[letter].dup }
    return PM.new(@size, new_matrix, @words_count)
  end
  
  def PM.new_pcm(words, iupacomp = false)
    size = words[0].size
    counts = PM.new_matrix(size)
    counts.each_value { |arr| arr.fill(0) }
    words.each { |word|
      0.upto(size-1) { |i|
        letter = word[i,1].upcase
        checkerr("unknown letter #{letter}") { !['A', 'C', 'G', 'T', 'N'].include?(letter) }
        if letter != 'N'
          counts[letter][i] += 1
        else
          ['A', 'C', 'G', 'T'].each { |l| counts[l][i] += 0.25 }
        end
      }
    }
    newpcm = PM.new(size, counts, words.size)
    newpcm.iupacomp! if iupacomp
    return newpcm
  end
  
  def PM.new_pwm(words)
    pcm = PM.new_pcm(words)
    pcm.to_pwm!
    return pcm
  end
  
  def PM.load(filename)
    # supporting pat & pwm formats (letter-column and letter-row format)
    input = IO.read(filename)
    tm = []
    input.each_line { |line|
      next if line.strip.empty?
      l_a = line.split
      begin
        l_a = l_a.collect { |a_i| Float(a_i) }
      rescue
        next
      end
      tm << l_a
    }
    tm = tm.transpose if tm.size == 4
    matrix = PM.new_matrix(tm.size)
    tm.each_index { |i| ['A', 'C', 'G', 'T'].each_with_index { |l, j| matrix[l][i] = tm[i][j] }  }
    
    ppm_mode = (0...tm.size).inject(true) { |ppm_ya, i| ppm_ya &= col_sum(matrix, i).round == 1 }
    
    return ppm_mode ? PPM.new(tm.size, matrix) : PM.new(tm.size, matrix)
  end
  
  def save(filename)
    File.open(filename, "w") { |out_f|
      case File.ext_wo_name(filename)
      when "pwm","ppm","pcm"
        ['A', 'C', 'G', 'T'].each { |letter|
          @matrix[letter].each { |e|
            out_f << "#{e} "
          }
          out_f << $/
        }
      when "pat"
        out_f.puts File.name_wo_ext(filename)
        (0...@size).each { |i|
          ['A', 'C', 'G', 'T'].each { |letter|
            out_f << "#{@matrix[letter][i]} "
          }
          out_f << $/
        }
      when "xml"
        checkerr("small-BiSMark is not supported at this moment")
      else
        checkerr("unknown motif file format specified")
      end
    }
  end
  
  def positiv!
    min = @matrix.values.collect { |v| v.min }.min.abs
    @matrix.each_value { |v| (0...v.size).each { |i| v[i] += min } }
    return self
  end
  
  def revcomp!
    @matrix['A'], @matrix['T'] = @matrix['T'], @matrix['A']
    @matrix['C'], @matrix['G'] = @matrix['G'], @matrix['C']
    @matrix.each_value { |v| v.reverse! }
    self
  end
  
  def to_bismark(b)
    pwm = @matrix['A'][0].is_a?(Float)
    attributes = {"length" => @size}
    attributes["words-count"] = @words_count if @words_count && @words_count > 0
    pe = b.add_element( pwm ? "PWM" : "PCM", attributes )
    (0...@matrix['A'].size).each { |i|
      pm_c = pe.add_element("pm-column", {"position" => i+1})
      ['A', 'C', 'G', 'T'].each { |l|
        pm_c.add_element(l.downcase).add_text(@matrix[l][i].to_s)
      }
    }
  end
  
  def PM.from_bismark(b, iupacomp = false)
    
    checkerr("empty small-BiSMark file?") { !b }
    float_m = (b.name == "PPM" || b.name == "PWM" || b.name == "WPCM")
    words_count = b.attributes["words-count"] ? b.attributes["words-count"].to_f : nil
    
    matrix = {"A" => [], "C" => [], "G" => [], "T" => []}
    b.elements.each("pm-column") { |pmc|
      position = pmc.attributes["position"].to_i
      ['A', 'C', 'G', 'T'].each { |l| matrix[l][position-1] = float_m ? pmc.elements[l.downcase].get_text.to_s.to_f : pmc.elements[l.downcase].get_text.to_s.to_i }
    }
    if b.name == "PPM"
      newppm = PPM.new(matrix['A'].size, matrix, words_count)
      newppm.iupacomp! if iupacomp
      return newppm
    end
    if b.name == "PCM"
      @words_count = col_sum(matrix)
      newpcm = PM.new(matrix['A'].size, matrix, words_count)
      newpcm.iupacomp! if iupacomp
      return newpcm
    end
    if b.name == "PWM" && iupacomp
      raise "cannot force IUPAC compatible PWM"
    end
    return PM.new(matrix['A'].size, matrix, words_count)
  end
  
  IUPAC_LS = (IUPAC::CODE.keys - ['A','C','G','T']).collect { |iul| [IUPAC::CODE[iul], iul.split(//)] }
  def iupacomp!
    @words_count = (0...@size).collect { |i| col_sum(i) }.max unless @words_count # for unbalanced matrices (Genomatix has some)
    # @words_count = @words_count.round < 2.0 ? nil : @words_count.round
    
    IUPAC_LS.each { |iul_ls|
      @matrix[iul_ls[0]] = (0...@size).collect { |i| col_sum(i, iul_ls[1]) / iul_ls[1].size }
    }
    
    return self
  end
  
  def m3sd(bckgr = Randoom::DEF_PROBS)
  
    mean = (0...@size).inject(0.0) { |mean, i| mean += ['A','C','G','T'].inject(0.0) { |sum,l| sum += @matrix[l][i] * bckgr[l] } }
    dev = (0...@size).inject(0.0) { |m2, i| 
      deltai = ['A','C','G','T'].inject(0.0) { |sum,l| sum += @matrix[l][i]**2 * bckgr[l] } - ['A','C','G','T'].inject(0.0) { |sum,l| sum += matrix[l][i] * bckgr[l] }**2
      m2 += deltai
    }
    sigma = Math.sqrt(dev)
    
    mean+3*sigma
  end
  
  def fixwc
    return unless @words_count
    @words_count = (0...@size).collect { |i| col_sum(i) }.max
  end
  
  protected
  def PM.new_matrix(size)
    return {
      'A' => Array.new(size),
      'C' => Array.new(size),
      'G' => Array.new(size),
      'T' => Array.new(size) }
  end
  
  def PM.new_matrix_iupac(size)
    return {
      'A' => Array.new(size),
      'C' => Array.new(size),
      'G' => Array.new(size),
      'T' => Array.new(size),
      'R' => Array.new(size),
      'Y' => Array.new(size),
      'K' => Array.new(size),
      'M' => Array.new(size),
      'S' => Array.new(size),
      'W' => Array.new(size),
      'B' => Array.new(size),
      'D' => Array.new(size),
      'H' => Array.new(size),
      'V' => Array.new(size),
      'N' => Array.new(size)
      }
  end
  
end

class PPM < PM
  
  #DEPRECATED, use iupacomp! instead
  #def make_N_comp!
  #  @matrix['N'] = (0...size).collect { 0.25 }
  #  return self
  #end
  
  def initialize(size, matrix = nil, words_count = nil)
    checkerr("matrix['A'].size != size") { matrix != nil && size != matrix['A'].size }
    @size = size
    @matrix = matrix == nil ? PM.new_matrix(size) : matrix
    @words_count = words_count
  end
  
  def iupacomp!
    @words_count = 4.0 unless @words_count
    
    IUPAC_LS.each { |iul_ls|
      @matrix[iul_ls[0]] = (0...@size).collect { |i| col_sum(i, iul_ls[1]) / iul_ls[1].size }
    }
    
    return self
  end
  
  def score(word)
    checkerr("word size != ppm.size") { @size != word.size }
    checkerr("word #{word} has strange characters") { 
      @matrix.keys.include?('N') ? word.tr('ACGTRYKMSWBDHVN', '').size > 0 : word.tr('ACGT', '').size > 0
    }
    return (0...@size).inject(1) { |mul, i| 
      mul *= @matrix[word[i,1]][i]
    }
  end
  
  def best_score
    return (0...size).inject(1) { |mul, i|
      mul *= ['A', 'C', 'G', 'T'].collect { |l| @matrix[l][i] }.max
    }    
  end
	
  def worst_score
    return (0...size).inject(0) { |mul, i|
      mul *= ['A', 'C', 'G', 'T'].collect { |l| @matrix[l][i] }.min
    }
  end
  
  def to_bismark(b)
    attributes = {"length" => @size}
    attributes["words-count"] = @words_count if @words_count
    pe = b.add_element("PPM", attributes)
    (0...@matrix['A'].size).each { |i|
      pm_c = pe.add_element("pm-column", {"position" => i+1})
      ['A', 'C', 'G', 'T'].each { |l|
        pm_c.add_element(l.downcase).add_text(@matrix[l][i].to_s)
      }
    }
  end
  
  def PPM.probs2IUPAC!(probs)
    IUPAC_LS.each { |iul_ls|
      probs[iul_ls[0]] = iul_ls[1].inject(0) { |sum, l| sum += probs[l] } / iul_ls[1].size
    }
    return probs
  end
  
  def get_pwm(words_count = nil, probs = Randoom::DEF_PROBS, pseudocount = 1.0)
    
    probs = PPM.probs2IUPAC!(probs.dup)
    
    words_count = @words_count if !words_count || words_count == 0
    checkerr("undefined words count") { !words_count }
    
    pwm = @matrix['N'] ? PM.new_matrix_iupac(@size) : PM.new_matrix(@size)
    
    @matrix.each_key do |letter|
      (0...@size).each { |pos|
        
        pwm[letter][pos] = Math::log( (@matrix[letter][pos] * words_count + (probs[letter] * pseudocount) ) / ( (words_count + pseudocount) * probs[letter]) )
        
      }
    end
    return PM.new(@size, pwm, words_count)
    #pcm = get_pcm(words_count)
    #pcm.iupacomp! if @matrix['N']
    #return pcm.to_pwm!(words_count, probs, pseudocount)
  end
  alias to_pwm get_pwm
  
  def get_pwm0pc(probs = Randoom::DEF_PROBS)
    new_matrix = {}
    @matrix.each_key { |letter| new_matrix[letter] = @matrix[letter].dup }
    newpm = PM.new(@size, new_matrix, nil)
    
    new_matrix.each_key do |letter|
      (0...@size).each { |pos|
        new_matrix[letter][pos] = Math::log(@matrix[letter][pos] / probs[letter])
      }
    end
    
    return newpm
  end
  
  def to_pwm!
    raise "cannot force PPM class to PWM, use to_pwm instead"
  end
  
  def get_pcm(words_count = nil)
    words_count = @words_count unless words_count
    checkerr("undefined words count") { !words_count }
    counts = PM.new_matrix(@size)
    (0...size).each { |i|
      ['A', 'C', 'G', 'T'].each { |l|
        counts[l][i] = @matrix[l][i] * words_count
      }
    }
    newpcm = PM.new(size, counts, words_count).iupacomp!
    pm = newpcm
    
    #(0...pm.size).each { |i|
    #  p pm.matrix['A'][i] + pm.matrix['C'][i] + pm.matrix['G'][i] + pm.matrix['T'][i]
    #}
    #p @words_count
    #p newpcm.words_count

    return newpcm
  end
  alias to_pcm get_pcm
  
  def PPM.from_IUPAC(iupac)
    matrix = {"A" => [], "C" => [], "G" => [], "T" => []}
    
    (0...iupac.size).each { |i|
      matrix.each_key { |k| matrix[k] << 0.0 }
      letters = IUPAC::REVCODE[iupac[i]]
      (0...letters.size).each { |j|
        matrix[letters[j]][-1] = 1.0/letters.size
      }
    }
    
    newppm = PPM.new(iupac.size, matrix, 4.0)
    newppm.iupacomp!
    
    newppm
  end
  
end

end
