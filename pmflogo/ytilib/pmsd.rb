#!/usr/bin/ruby
module Ytilib
  class PM
    def score_sigma(trycount = 4**10, approx = false, bg = nil)
      
      scores = []
      if @size <= 10 && !approx
        (0...4**@size).each { |i| 
          word = i.to_s(4).rjust(@size, "0").tr("0123", "ACGT")
          scores << score(word)
        }
      else
        trycount.times {
          word = bg ? Randoom.rand_seq(@size, bg) : Randoom.rand_seq(@size)
          scores << score(word)
        }
      end
      sum1 = scores.inject(0) { |sum,s| sum += s }
      mean = sum1 / scores.size
      
      sum2, sumc = 0, 0
      scores.each { |score|
        sum2 += (score-mean)**2
        sumc += (score-mean)
      }
      variance = (sum2 - sumc**2 / scores.size) / (scores.size-1)
      
      sigma = Math.sqrt(variance)
      if block_given?
        yield(sigma, mean)
      end
      
      return sigma
    end
    
    def fast_score_sigma
      n, mean, m2 = 0, 0, 0
      
      recursive_walk([matrix['A'],matrix['C'],matrix['G'],matrix['T']], 0, 0) { |x|
        n = n + 1
        delta = x - mean
        mean = mean + delta/n
        m2 = m2 + delta*(x-mean)
      }
      
      variance = m2/(n - 1)
      
      if block_given?
        yield(sigma = Math.sqrt(variance), mean)
      end
      
      return sigma
    end
    
    def fast_score_sigma_precise
      n, mean = 0, 0
      
      recursive_walk([matrix['A'],matrix['C'],matrix['G'],matrix['T']], 0, 0) { |x|
        n += 1
        delta = x - mean
        mean = mean + delta/n
      }
      
      n, m2 = 0, 0
      recursive_walk([matrix['A'],matrix['C'],matrix['G'],matrix['T']], 0, 0) { |x|
        n = n + 1
        delta = x - mean
        m2 = m2 + delta*(x-mean)
      }
      
      variance = m2/(n - 1)
      
      if block_given?
        yield(sigma = Math.sqrt(variance), mean)
      end
      
      return sigma
    end
    
  private
    def recursive_walk(matrix, score, i)
      if i < @size
        
        recursive_walk(matrix, score + matrix[0][i], i+1) { |x| yield x }
        recursive_walk(matrix, score + matrix[1][i], i+1) { |x| yield x }
        recursive_walk(matrix, score + matrix[2][i], i+1) { |x| yield x }
        recursive_walk(matrix, score + matrix[3][i], i+1) { |x| yield x }
        
      else
        if block_given?
          yield(score)
        else
          raise "no block for recursive walk"
        end
      end
    end
  
  end
end