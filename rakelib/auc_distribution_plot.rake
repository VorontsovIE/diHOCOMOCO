task :auc_distribution_plot do
  auc_by_control_and_model = Hash.new{|h,k| h[k] = {} }
  Models.all_models.each{|model|
    model_auc_fn = File.join('occurences/auc/', model.uniprot, "#{model.full_name}.txt")
    File.readlines(model_auc_fn).each{|line|
      control_name, auc = line.chomp.split("\t")
      auc_by_control_and_model[control_name][model.full_name] = auc.to_f
    }
  }

  max_aucs_for_control = auc_by_control_and_model.map{|control, auc_by_model|
    auc_by_model.values.max
  }

  auc_distribution = (0.0..1.0).step(0.01).map{|t|
    cnt = max_aucs_for_control.count{|x| x <= t }
    [t, cnt]
  }
  File.write 'auc_distribution.tsv', auc_distribution.map{|t, cnt| "#{t}\t#{cnt}" }.join("\n")
end
