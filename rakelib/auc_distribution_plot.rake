require 'auc_by_model_and_control'
task :auc_distribution_plot do
  max_aucs_for_control = auc_by_control_and_model('occurences/auc/').map{|control, auc_by_model|
    auc_by_model.values.max
  }

  auc_distribution = (0.0..1.0).step(0.01).map{|t|
    cnt = max_aucs_for_control.count{|x| x <= t }
    [t, cnt]
  }
  File.write 'auc_distribution.tsv', auc_distribution.map{|t, cnt| "#{t}\t#{cnt}" }.join("\n")
end
