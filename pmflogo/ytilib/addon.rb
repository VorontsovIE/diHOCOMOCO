#!/usr/bin/ruby

def File.ext_wo_name(what)
  return what if what.rindex(".") == nil
  what = File.basename(what)
  "#{what}"[what.rindex(".")+1..-1]
end

def File.name_wo_ext(what)
  return what if what.rindex(".") == nil
  what = File.basename(what)
  "#{what}"[0...what.rindex(".")]
end

class Float
  def round_to(x)
    (self * 10**x).round.to_f / 10**x
  end
  
  def cut_to(x)
    (self.abs * 10**x).floor.to_f * (self == 0.0 ? 0 : self/self.abs).round / 10**x
  end
end

class Array
  def shuffle
    arr = self.dup
    arr.size.downto 2 do |j|
      r = rand j
      arr[j-1], arr[r] = arr[r], arr[j-1]
    end
    arr
  end
  
  def shuffle!
    (size - 1).downto 1 do |i|
      j = rand(i + 1)
      self[i], self[j] = self[j], self[i]
    end
    self
  end
  
  def average
    self.empty? ? nil : self.inject(0) { |sum,s| sum += s } / self.size
  end
  alias mean average
  
  def variance
    return self.collect { |s| s*s }.average - average**2
  end
  
  def sum
    self.inject(self[0]) { |sum,s| sum += s} - self[0]
  end

end

class String
  
  def compl!
    self.tr!("acgtACGT", "tgcaTGCA")
    return self
  end
  
  def compl
    return self.tr("acgtACGT", "tgcaTGCA")
  end
  
  alias comp! compl!
  alias complement! compl!
  alias comp compl
  alias complement compl
  
  def revcomp
    return comp.reverse
  end
  
  def revcomp!
    return comp!.reverse!
  end
  
  def to_id
    return self.gsub(/[^.\w]/, '_').upcase
  end
  
end

# Also this can be done is a more sophisticated way
=begin
String.class_eval do    
  def to_id
    return self.gsub(/[^.\w]/, '_') 
  end   
end
=end

class String
  # The opposite of String::next / String::succ. It is impossible to be a
  # *complete* opposite because both "9".next = "10" and "09".next = "10";
  # if going backwards from "10" there's no way to know whether the result
  # should be "09" or "9". Where the first ranged character is about to
  # underflow and the next character is within the same range the result
  # is shrunk down - that is, "10" goes to "9", "aa" goes to "z"; any non-
  # range prefix or suffix is OK, e.g. "+!$%10-=+" goes to "+!$%9-=+".
  # Items in the middle of a string don't do this - e.g. "12.10" goes to
  # "12.09", to match how "next" would work as best as possible.
  #
  # The standard "next" function works on strings that contain *no*
  # alphanumeric characters, using character codes. This implementation
  # of "prev" does *not* work on such strings - while strings may contain
  # any characters you like, only the alphanumeric components are operated
  # upon.
  #
  # Should total underflow result, "nil" will be returned - e.g. "00".prev
  # returns 'nil', as does "a".prev. This is done even if there are other
  # characters in the string that were not touched - e.g. "+0.0".prev
  # also returns "nil". Broadly speaking, a "nil" return value is used for
  # any attempt to find the previous value of a string that could not have
  # been generated using "next" in the first place.
  #
  # As with "next" sometimes the result of "prev" can be a little obscure
  # so it is often best to try out things using "irb" if unsure. Note in
  # particular that software revision numbers do not necessarily behave
  # predictably, because they don't with "next"! E.g. "12.4.9" might go to
  # "12.4.10" for a revision number, but "12.4.9".next = "12.5.0". Thus
  # "12.5.0".prev = "12.4.9" and "12.4.10".prev = "12.4.09" (because the
  # only way to make "12.4.10" using "next" is to start at "12.4.09").
  #
  # Since 'succ' (successor) is an alias for 'next', so 'pred'
  # (predecessor) is an alias for 'prev'.
  #
  def prev(collapse = false)
    str        = self.dup
    early_exit = false
    any_done   = false
    ranges     = [
                   ('0'[0]..'9'[0]),
                   ('a'[0]..'z'[0]),
                   ('A'[0]..'Z'[0]),
                   nil
                 ]

    # Search forward for the first in-range character. If found check
    # to see if that character is "1", "a" or "A". If it is, record
    # its index (from 0 to string length - 1). We'll need this if
    # underflows wrap as far as the found byte because in that case
    # this first found byte should be deleted ("aa..." -> "z...",
    # "10..." -> "9...").

    first_ranged = nil

    for index in (1..str.length)
      byte = str[index - 1]

      # Determine whether or not the current byte is a number, lower case
      # or upper case letter. We expect 'select' to only find one matching
      # array entry in 'ranges', thus we dereference index 0 after the
      # 'end' to put a matching range from within 'ranges' into 'within',
      # or 'nil' for any unmatched byte.

      within = ranges.select do |range|
        range.nil? or range.include?(byte)
      end [0]

      unless within.nil?
        case within.first
          when '0'[0]
            match_byte = '1'[0]
          else
            match_byte = within.first
        end

        first_ranged = index - 1 if (byte == match_byte)
        first_within = within
        break
      end
    end

    for index in (1..str.length)

      # Process the input string in reverse character order - fetch the
      # bytes via negative index.

      byte = str[-index]

      within = ranges.select do |range|
        range.nil? or range.include?(byte)
      end [0]

      # Skip this letter unless within a known range. Otherwise note that
      # at least one byte was able to be processed.

      next if within.nil?
      any_done = true

      # Decrement the current byte. If it is still within its range, set
      # the byte and bail out - we're finished. Flag the early exit. If
      # the byte is no longer within range, wrap the character around
      # and continue the loop to carry the decrement to an earlier byte.

      byte = byte - 1

      if (within.include? byte)
        str[-index] = byte
        early_exit  = true
        break
      else
        str[-index] = within.last

        # If we've just wrapped around a character immediately after the
        # one found right at the start ('0', 'a' or 'A') then this first
        # ranged character should be deleted (so "10" -> "09"

        if (first_ranged != nil and first_within.include?(byte + 1) and (first_ranged - str.length) == -(index + 1))
          str.slice!(-(index + 1))
          early_exit = true
          break
        end
      end

    end # From outer 'for' loop

    # If we did process at least one byte but we did not exit early, then
    # the loop completed due to carrying a decrement to other bytes. This
    # means an underflow condition - return 'nil'.

    if (any_done == true and early_exit == false)
      return nil
    else
      return str
    end
  end

  # As (extended) String::pred / String::prev, but modifies the string in
  # place rather than returning a copy. If underflow occurs, the string
  # will be unchanged. Returns 'self'.
  #
  def prev!
    new_str = prev
    self.replace(new_str) unless new_str.nil?
    return self
  end

  alias pred  prev
  alias pred! prev!

end