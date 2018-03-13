['mono', 'di'].each do |model_type|
  task "wlogaucs_for_slices_#{model_type}" do
    semiuniprots = Dir.glob("auc/#{model_type}/*_datasets/*").map{|fn|
      File.basename(fn)
    }.map{|basename|
      basename.split('~').first.split('_').first
    }.uniq.sort

    semiuniprots.each do |semiuniprot|
      puts 'ruby weight_aucs.rb #{model_type} #{semiuniprot}'
    end
  end
end
