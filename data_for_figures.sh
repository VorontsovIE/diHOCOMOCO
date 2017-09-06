mkdir -p figures
find control/control/ -xtype f | ruby -e 'result = readlines.map{|fn| File.basename(fn.chomp, ".control.mfa") }.group_by{|peakset| peakset.split(".")[0] }.map{|tf,peaksets| datasets = peaksets.map{|peakset| peakset.split(".")[1] }.uniq; [tf, peaksets.size, datasets.size] }; result.each{|infos| puts infos.join("\t") }' > figures/num_peaksets_and_datasets_by_tf.tsv
