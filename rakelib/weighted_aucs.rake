require 'fileutils'

['mono', 'di'].each do |model_type|
  task "wlogaucs_for_slices_#{model_type}" do
    FileUtils.rm_r "wlogauc/#{model_type}"  if Dir.exist?("wlogauc/#{model_type}")

    semiuniprots = Dir.glob("auc/#{model_type}/*_datasets/*").map{|fn|
      File.basename(fn)
    }.map{|basename|
      basename.split('~').first.split('_').first
    }.uniq.sort

    semiuniprots.each do |semiuniprot|
      puts "ruby weight_aucs.rb #{model_type} #{semiuniprot}"
    end
  end
end

task "calculate_wlogaucs" => ["wlogaucs_for_slices_mono", "wlogaucs_for_slices_di"]
