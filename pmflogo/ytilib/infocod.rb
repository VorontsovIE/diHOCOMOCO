#!/usr/bin/ruby

class Float
  # Using Stieltjes formula from http://www.luschny.de/math/factorial/approx/SimpleCases.html
  def log_fact
    return 0.0 if self <= 1
    a0 = 1.0/12
    a1 = 1.0/30 
    a2 = 53.0/210
    a3 = 195.0/371
    a4 = 22999.0/22737
    a5 = 29944523.0/19733142
    a6 = 109535241009.0/48264275462
    z_big = self+1;
    (1.0/2)*Math.log(2*Math::PI)+(z_big-1.0/2)*Math.log(z_big)-z_big + a0/(z_big+a1/(z_big+a2/(z_big+a3/(z_big+a4/(z_big+a5/(z_big+a6/z_big))))))
  end
end

class Integer
  def log_fact
    self.to_f.log_fact
  end
end

# Naive version
=begin
class Integer
  @@fact_hash = {}
  def log_fact
    return 0.0 if self == 0
    return nil if self < 0
    if self <= 170
      @@fact_hash[self] = Math.log( lambda { |k| return k if self.times { |i| k *= i.next } }.call(1) )
    else
      return self.to_f.log_fact
    end unless @@fact_hash.has_key?(self)
    return @@fact_hash[self] 
  end
end
=end

module Ytilib
  class PM
    def infocod(position = nil)
      return infocod_private(position) if position
      (0...@size).collect { |i| infocod_private(i) }
    end
    
    def icd2of4(floor = false)
      i2o4 = @words_count / 2.0
      i2o4 = i2o4.floor if floor
      ([i2o4, i2o4, 0, 0].inject(0.0) { |sum, k_i| sum += k_i.log_fact  } - @words_count.log_fact ) / @words_count
      # 0 is equal to @words_count % 2, because 0! = 1!
    end
    
    def icd3of4(floor = false)
      i3o4 = @words_count / 3.0
      i3o4 = i3o4.floor if floor
      addon = floor ? @words_count % 3 : 0
      ([i3o4, i3o4, i3o4, addon].inject(0.0) { |sum, k_i| sum += k_i.log_fact  } - @words_count.log_fact ) / @words_count
    end
    
    def icdThc
      icd3of4
    end
    
    def icdTlc
      io = @words_count / 6.0
      ([2*io, 2*io, io, io].inject(0.0) { |sum, k_i| sum += k_i.log_fact  } - @words_count.log_fact ) / @words_count
    end
    
    def icd4of4(floor = false)
      i4o4 = @words_count / 4.0
      i4o4 = i4o4.floor if floor
      ([i4o4, i4o4, i4o4, i4o4].inject(0.0) { |sum, k_i| sum += k_i.log_fact  } - @words_count.log_fact ) / @words_count
    end
    
  protected
    def infocod_private(position)
      k_i = ['A','C','G','T'].collect { |letter| @matrix[letter][position] }
      
      ( k_i.inject(0.0) { |sum, k_i| sum += k_i.log_fact  } - @words_count.log_fact ) / @words_count
    end
  end
  
  class PPM
    #def to_pcm(words_count = nil)
     # @words_count = words_count if words_count
     # checkerr("words count is not specified") { !@words_count }
     # counts = PM.new_matrix(@size)
     # (0...size).each { |i|
     #   ['A', 'C', 'G', 'T'].each { |l|
     #     counts[l][i] = @matrix[l][i] * @words_count
     #   }
     # }
     # return PM.new(size, counts)
     #end
     #alias to_pcm get_pcm
    
    def infocod(position = nil)
      return to_pcm.infocod(position)
    end
    alias dic infocod
    
    def icd(position = nil)
      return to_pcm.infocod(position)
    end
  end
end
