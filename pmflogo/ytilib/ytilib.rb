#!/usr/bin/ruby

module Ytilib
  def Ytilib.time
    return Time.now.strftime('%d %b %H:%M:%S')
  end
end

$program_name = nil

def start(fullpath)
  report(fullpath + ARGV.inject("") { |out, v| out += " " + v})
  return if $NO_REPORT
  $program_name = "[#{File.name_wo_ext(fullpath)}]"
end

def report(message, program_name = nil)
  $program_name = "[#{program_name}]" if program_name != nil
  return if $NO_REPORT
  puts "LLIB #{Ytilib.time} #{$program_name}\t#{message}" if !block_given? || yield
end

def checkerr(message = "checkerr failed")  
  if !block_given? || yield
    puts "LLIB #{Ytilib.time} [error]\t#{message}" unless $NO_REPORT
    raise "LLIB #{Ytilib.time} #{$program_name}\n\t#{message}\n" 
  end
end

module Ytilib
  
  STRAND_DIRECT = "direct"
  STRAND_REVCOMP = "revcomp"
  
  #require 'seqplace.rb'  
  #require 'mysql.rb'
  require 'addon.rb'
  require 'iupac.rb'
  require 'pm.rb'
  require 'pmsd.rb'
  require 'randoom.rb'
  require 'bismark.rb'
  require 'infocod.rb'
  require 'hack1.rb'
  
  def Ytilib.read_mfa2hash(path)
    input_fasta_f = File.new(path, "r")
    seqs, seq_name = {}, nil
    input_fasta_f.each_line { |line|
      if line[0,1] == ">"
        seq_name = line[1..-1].strip
        seq_name = yield seq_name if block_given?
        checkerr("multiple sequences with the same name=#{seq_name}") { seqs[seq_name] }
        seqs[seq_name] = ""
      elsif seq_name != nil
        seqs[seq_name] << line.strip
      end
    }
    input_fasta_f.close
    return seqs
  end
  
  def Ytilib.read_mfa2array(path)
    input_fasta_f = File.new(path, "r")
    seqs, seq_name = [], nil
    input_fasta_f.each_line { |line|
      if line[0,1] == ">"
        seq_name = line[1..-1].strip
        yield seq_name if block_given?
        seqs << ""
      elsif seq_name != nil
        seqs.last << line.strip
      end
    }
    input_fasta_f.close
    return seqs
  end
  
  def Ytilib.mfa2array(input)
    seqs, seq_name = [], nil
    input.each_line { |line|
      if line[0,1] == ">"
        seq_name = line[1..-1].strip
        seqs << ""
      elsif seq_name != nil
        seqs.last << line.strip
      end
    }
    return seqs
  end
  
  def Ytilib.read_plain2array(path)
    array = []
    File.open(path).each_line { |line|
      array << line.strip if !line.strip.empty?
    }
    return array
  end  
  
  def Ytilib.read_seqs2array(path)
    type = File.ext_wo_name(path)
    case type
    when "mfa", "fasta", "fa"
      return Ytilib.read_mfa2array(path)
    when "plain","txt"
      return Ytilib.read_plain2array(path)
    else
      checkerr("unknown sequences-file, ext=#{type}")
    end
  end
  
  def Ytilib.write_mfa(seqs, path, prefix = " ")
    if seqs.is_a?(Hash)
      out_fasta_f = File.new(path, "w+")
      seqs.each_key { |name|
        out_fasta_f << ">#{prefix}#{name}" << $/ << seqs[name] << $/
      }
      out_fasta_f.close
    else 
      out_fasta_f = File.new(path, "w+")
      seqs.each_with_index { |seq, i|
        out_fasta_f << ">#{prefix}#{i+1}" << $/ << seq << $/
      }
      out_fasta_f.close
    end
  end
  
  def get_consensus(seqs)
    report "consensus creating method should be checked, you are using unsafe code"
    return 'nil' if seqs.size == 0
    conslet = { 'A' => 'A', 'C' => 'C', 'G' => 'G', 'T' => 'T', 'U' => 'U',
            'AG' => 'R', 'CT' => 'Y', 'GT' => 'K', 'AC' => 'M', 'CG' => 'S', 'AT' => 'W',
            'CGT' => 'B', 'AGT' => 'D', 'ACT' => 'H', 'ACG' => 'V', 'ACGT' => 'N'
    }
    new_consensus, letters = '', []
    0.upto(seqs[0].size-1) { |i|
      seqs.each do |word|
        letters << word[i] if !letters.include?(word[i])
      end
      letters.sort!
      letters_string = ''
      letters.each do |letter| letters_string << letter end
      checkerr("cannot find consensus letter for a given letter set :#{}") { conslet[letters_string] == nil }
      new_consensus << conslet[letters_string]	
      letters.clear
    }
    return new_consensus
  end

  def Ytilib.new_mysql_conn(database)
    my = Mysql.new(MYSQL_HOST, MYSQL_USER, MYSQL_PASSWORD, database)
    checkerr("cannot connect to MySQL server") { my.query("select 1").fetch_row[0] != "1" }
    return my
  end  
  

end

report "ytilib required, working directory #{Dir.pwd}", "ytilib"
include Ytilib
