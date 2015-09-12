#!/usr/bin/ruby
app_dir = File.dirname(File.expand_path(__FILE__))
$: << "#{app_dir}/ytilib" 
$: << "#{app_dir}/../lib" 

require "ytilib.rb" 
basef = File.dirname(__FILE__)

require 'rmagick'
require 'bioinform'
require_relative '../lib/models'
require_relative '../lib/dipm'


ns = ['A','C','G','T']
$dins = dins = ns.collect { |n| ns.collect { |n1| "#{n}#{n1}" } }.flatten



def load_dipm(path, should_revcomp:)
  if should_revcomp
    ModelKind.get('di').read_pcm(path).revcomp.to_hash
  else
    ModelKind.get('di').read_pcm(path).to_hash
  end
end

def sum_din_beg(letter, dipm, position)
  $dins.select { |d| d[0,1] == letter }.inject(0.0) { |sum,s| sum += dipm[s][position] }
end

def sum_din_end(letter, dipm, position)
  $dins.select { |d| d[1,1] == letter }.inject(0.0) { |sum,s| sum += dipm[s][position] }
end

def make_pm(dipm)
  pm = {'A' => [], 'C' => [], 'G' => [], 'T' => []}
  size = dipm['AA'].size-1
  
  # ['A','C','G','T'].each { |let| pm[let] << sum_din_beg(let, dipm, 0) }
  
  
  (0...size).each { |pos|
    ['A','C','G','T'].each { |let| pm[let] << sum_din_end(let, dipm, pos) }
  }
  pm
end

def make_pma(dipm)
  pm = {'A' => [], 'C' => [], 'G' => [], 'T' => []}
  size = dipm['AA'].size
  
  
  (1...size).each { |pos|
    ['A','C','G','T'].each { |let| pm[let] << sum_din_beg(let, dipm, pos) }
  }
  pm
end

report "pmflogo3.rb started, usage: <in_file_name.dpcm/fa/mfa> <out_file_name> [<x_unit>=100] [<y_unit>=200] [<mono_scheme>=nucl_simpa] [<di_scheme>=di_simpa]"
start __FILE__
exit(2) if ARGV.size < 2 

should_revcomp = ARGV.delete('--revcomp')

PAR_input = ARGV[0]
PAR_output = ARGV[1]
x_unit = ARGV[2] ? ARGV[2].to_i : 100
y_unit = ARGV[3] ? ARGV[3].to_i : 200
mono_scheme = ARGV[4] ? ARGV[4] : "nucl_simpa"
di_scheme = ARGV[5] ? ARGV[5] : "dinucl_simpa"

dipm = load_dipm(PAR_input, should_revcomp: should_revcomp)
pm = make_pm(dipm)
#pm = make_pma(dipm)

N = (0...dipm['AA'].size).collect { |pos|
  dins.inject(0.0) { |sum,s| sum += dipm[s][pos] }
}.max

def dikdic(n, vs, pos, background)
  r = n.log_fact
  $dins.each { |din|
    r -= vs[din][pos].log_fact
    r += vs[din][pos] * Math.log(background[din])
  }
  return -r / n
end

def kdic(n, vs, pos)
  r = n.log_fact
  ['A','C','G','T'].each { |nu|
    r -= vs[nu][pos].log_fact
    r += vs[nu][pos] * Math.log(0.25)
  }
  return -r / n
end

background = {}
dins.each { |d| background[d] = 0.0625 }

dikdic_r = (0...dipm['AA'].size).collect { |pos|
  dikdic(N, dipm, pos, background)
}

kdic_r = (0...pm['A'].size).collect { |pos|
  kdic(N, pm, pos)
}

def dikdic_max(n)
  #r = n.log_fact
  #r -= n.log_fact
  r = n * Math.log(0.0625)
  -r / n
end

def dikdic_min(n)
  r = n.log_fact
  $dins.each { |din|
    r -= (n/16.0).log_fact
    r += (n/16.0) * Math.log(0.0625)
  }
  return -r / n
end

def kdic_max(n)
  r = n*Math.log(0.25)
  -r / n
end

def kdic_min(n)
  r = n.log_fact
  ['A','C','G','T'].each { |nu|
    r -= (n/4.0).log_fact
    r += (n/4.0) * Math.log(0.25)
  }
  return -r / n
end

def dikdic_scale(v, min = dikdic_min(N), max = dikdic_max(N))
  (v-min) / (max-min)
end

def kdic_scale(v, min = kdic_min(N), max = kdic_max(N))
  (v-min) / (max-min)
end

def dikdic_min(n)
  r = n.log_fact
  $dins.each { |din|
    r -= (n/16.0).log_fact
    r += (n/16.0) * Math.log(0.0625)
  }
  return -r / n
end

report "KDIDIC = #{dikdic_r.inject(0) { |sum,s| sum += s } / (dikdic_r.size * dikdic_max(N))}"
report "weight = #{N}"

# prepare scale - column by dikdic, letters by counts (div by N)
dipm_scaled = {}
(0...dipm['AA'].size).each { |pos|
  dins.each { |din|
    dipm_scaled[din] = [] unless dipm_scaled[din]
    dipm_scaled[din][pos] = dikdic_scale(dikdic_r[pos]) * dipm[din][pos] / N
  }
}
pm_scaled = {}
(0...pm['A'].size).each { |pos|
  ns.each { |nu|
    pm_scaled[nu] = [] unless pm_scaled[nu]
    pm_scaled[nu][pos] = kdic_scale(kdic_r[pos]) * pm[nu][pos] / N
  }
}

#x_unit, y_unit = 100, 200
x_size, y_size = x_unit * dipm['AA'].size, 2*y_unit

i_dins = Magick::ImageList.new
dins.collect { |din| "#{basef}/#{di_scheme}/#{din.downcase}.png" }.each { |f| i_dins.read(f) }

i_ns = Magick::ImageList.new
ns.collect { |let| "#{basef}/#{mono_scheme}/#{let.downcase}.png" }.each { |f| i_ns.read(f) }

i_logo = Magick::ImageList.new
i_logo.new_image(x_size, y_size, Magick::HatchFill.new('white', 'white'))

letter_indexes = {}
dins.each_with_index { |d,i| letter_indexes[d] = i }
ns.each_with_index { |nu,i| letter_indexes[nu] = i }

(0...dipm['AA'].size).each { |i|
  
  #y_pos = 0
  y_pos = y_size/2.0
  sorted_dins = dins.collect { |letter| {:score => dipm[letter][i], :letter => letter} }.sort_by { |pair| pair[:score] }.collect { |pair| pair[:letter] }.reverse  
  
  sorted_dins.each { |letter|
    next if y_unit * dipm_scaled[letter][i] <= 1
    letter_index = letter_indexes[letter]
    y_block = (y_unit * dipm_scaled[letter][i]).round
    i_logo << i_dins[letter_index].dup.resize(x_unit, y_block)
    #i_logo.cur_image.page = Magick::Rectangle.new(0, 0, i*x_unit, y_size - y_size / 2.0 - y_pos )
    i_logo.cur_image.page = Magick::Rectangle.new(0, 0, i*x_unit, y_pos + 2)
    y_pos += y_block
  }
  
  #y_pos = y_size/2.0
  y_pos = 0
  sorted_ns = ns.collect { |letter| {:score => pm[letter][i], :letter => letter} }.sort_by { |pair| pair[:score] }.collect { |pair| pair[:letter] }.reverse unless i == pm['A'].size
  
  sorted_ns.each { |letter|
    next if y_unit * pm_scaled[letter][i] <= 1
    letter_index = letter_indexes[letter]
    y_block = (y_unit * pm_scaled[letter][i]).round
    i_logo << i_ns[letter_index].dup.resize(x_unit, y_block)
    #i_logo.cur_image.page = Magick::Rectangle.new(0, 0, i*x_unit + x_unit / 2.0, y_pos )
    
    y_pos += y_block
    i_logo.cur_image.page = Magick::Rectangle.new(0, 0, i*x_unit + x_unit / 2.0, y_size - y_size / 2.0 - y_pos - 1)
    
  } unless i == pm['A'].size
  
}

i_logo = i_logo.flatten_images
dr = Magick::Draw.new
dr.fill('transparent')
dr.stroke_width(y_size / 300.0)
dr.stroke_dasharray(20,3)

#dr.stroke('lightblue')
dr.stroke('#d9d9d9')
#dr.stroke('gray')
dr.line(0, y_size/2.0, x_size, y_size/2.0)

(0...dipm['AA'].size).each { |i|
  dr.line(i*x_unit + x_unit / 2.0, 0, i*x_unit + x_unit / 2.0, y_size / 2.0)
}
(0...dipm['AA'].size).each { |i|
  dr.line(i*x_unit, y_size / 2.0, i*x_unit, y_size)
}

dr.draw(i_logo)
  
i_logo.write(PAR_output)
