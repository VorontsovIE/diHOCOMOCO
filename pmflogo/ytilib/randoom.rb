#!/usr/bin/ruby
module Ytilib

srand

module Randoom
  
  private

  def Randoom.new_counts
    { 'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0, 'N' => 0 }
  end
  
  def Randoom.random_letter(probs)
    random = rand()
    return 'A' if random < probs['A']
    return 'C' if random < probs['A'] + probs['C']
    return 'G' if random < probs['A'] + probs['C'] + probs['G']
    return 'T'
  end

  public
  
  def Randoom.calc_probs(input)
    counts = new_counts
    counts.default = 0
    (0...input.length).each { |i|
      counts[input[i,1].upcase] += 1
    }
    return make_probs!(counts)
  end
    
  def Randoom.rand_seq(req_len, probs = DEF_PROBS, probs_m = nil)
    randoom = ''
    if (probs_m == nil)
      req_len.times { randoom << random_letter(probs) }
      return randoom  
    end
    random_l = random_letter(probs)
    randoom = random_l
    (req_len-1).times {
      cur_probs = probs_m[random_l]
      random_l = random_letter(cur_probs)
      randoom << random_l
    }
    return randoom
  end
  
  def Randoom.calc_probs_m(input)
    probs_m = { 'A' => {}, 'C' => {}, 'G' => {}, 'T' => {} }
    counts = { 'A' => new_counts, 'C' => new_counts, 'G' => new_counts, 'T' => new_counts, 'N' => new_counts }
    (0...input.length-1).each { |i|
      pair = input[i, 2].upcase
      counts[pair[0,1]][pair[1,1]] += 1
    }
    probs_m['A'] = make_probs!(counts['A'])
    probs_m['C'] = make_probs!(counts['C'])
    probs_m['G'] = make_probs!(counts['G'])
    probs_m['T'] = make_probs!(counts['T'])
    return probs_m
  end
  
  def Randoom.make_probs_m!(counts)
    ['A','C','G','T','N'].each { |l2|
      addv = counts['N'][l2] / 4.0
      ['A','C','G','T'].each { |l1|
        counts[l1][l2] += addv
      }
    }
    
    probs_m = { 'A' => {}, 'C' => {}, 'G' => {}, 'T' => {} }
    probs_m['A'] = make_probs!(counts['A'])
    probs_m['C'] = make_probs!(counts['C'])
    probs_m['G'] = make_probs!(counts['G'])
    probs_m['T'] = make_probs!(counts['T'])
    return probs_m
  end
  
  def Randoom.make_probs!(counts, length = nil)
    probs = { 'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0 }
    length = counts['A'] + counts['C'] + counts['G'] + counts['T'] + counts['N'] if length == nil
    length = length.to_f
    ['A','C','G','T'].each { |l| counts[l] += counts['N'] / 4.0 }
    return probs if length == 0
    probs['A'] = counts['A'] / length
    probs['C'] = counts['C'] / length
    probs['G'] = counts['G'] / length
    probs['T'] = 1 - probs['A'] - probs['C'] - probs['G']
    return probs
  end
  
  def Randoom.equalize!(probs)
    probs['A'] = probs['T'] = (probs['A'] + probs['T']) / 2
    probs['C'] = probs['G'] = (probs['C'] + probs['G']) / 2
    return probs
  end
  
  def Randoom.twostrand!(probs)
    return Randoom.equalize!(probs)
  end
  
  DEF_PROBS = PPM.probs2IUPAC!({ 'A' => 0.25, 'C' => 0.25, 'G' => 0.25, 'T' => 0.25, 'N' => 0.25 })
  
  # probabilities counted without _random.fa files for human genome
  DMEL40_PROBS1 = {"A"=>0.287729562173578, "C"=>0.21236364146414, "G"=>0.212259972960341, "T"=>0.287646823401942}
  DMEL40_PROBS2 = {"A"=>0.28768819278776, "C"=>0.21231180721224, "G"=>0.21231180721224, "T"=>0.28768819278776}
  
  DMEL40_PROBS1_M = {"A"=>{"A"=>0.350403075314602, "C"=>0.181194374386404, "G"=>0.188361404205017, "T"=>0.280041146093977}, 
    "C"=>{"A"=>0.325366772443085, "C"=>0.222264645612127, "G"=>0.197213801868993, "T"=>0.255154780075794}, 
    "G"=>{"A"=>0.260710563672393, "C"=>0.27150575901391, "G"=>0.222294234776053, "T"=>0.245489442537644}, 
    "T"=>{"A"=>0.217189093089999, "C"=>0.192590127484359, "G"=>0.239869076706963, "T"=>0.350351702718679}}
  
  HG17_PROBS1 = {"A"=>0.295309361730334, "C"=>0.204413561169847, "G"=>0.204519414193999, "T"=>0.295757662905821}
  HG17_PROBS2 = {"A"=>0.295533512318077, "C"=>0.204466487681923, "G"=>0.204466487681923, "T"=>0.295533512318077}
  
  HG17_PROBS1_M = {"A"=>{"A"=>0.331091206257755, "C"=>0.170458424092748, "G"=>0.236770972081246, "T"=>0.261679397568252}, 
    "C"=>{"A"=>0.354813019140533, "C"=>0.254741288394943, "G"=>0.0481667110625576, "T"=>0.342278981401966}, 
    "G"=>{"A"=>0.290057117684408, "C"=>0.208514091370804, "G"=>0.254732297362797, "T"=>0.246696493581991}, 
    "T"=>{"A"=>0.222087715262152, "C"=>0.200697606508443, "G"=>0.245657322003887, "T"=>0.331557356225517}}
  
  HG18_PROBS1 = {"A"=>0.291900580635872, "C"=>0.207855064518284, "G"=>0.207968587245859, "T"=>0.292275767599985}
  HG18_PROBS2 = {"A"=>0.292088174117929, "C"=>0.207911825882071, "G"=>0.207911825882071, "T"=>0.292088174117929}
  
  MM9_PROBS1 = {"A"=>0.289755259854654, "C"=>0.210085673636132, "G"=>0.210143929198141, "T"=>0.290015137311074}
  MM9_PROBS2 = {"A"=>0.289885198582864, "C"=>0.210114801417136, "G"=>0.210114801417136, "T"=>0.289885198582864}
  
  MM9_PROBS1_M = {"A"=>{"A"=>0.310389104265713, "C"=>0.184962574392377, "G"=>0.251904718465914, "T"=>0.252743602875996}, "C"=>{"A"=>0.352189584318682, "C"=>0.250794045222924, "G"=>0.0494404816637487, "T"=>0.347575888794645}, "G"=>{"A"=>0.295931117515178, "C"=>0.197870954111653, "G"=>0.250756985626016, "T"=>0.255440942747154}, "T"=>{"A"=>0.219437756702452, "C"=>0.214548041970626, "G"=>0.255405334730743, "T"=>0.310608866596179}}

end

end
