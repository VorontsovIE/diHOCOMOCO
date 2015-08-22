require 'auc_by_model_and_control'
task :auc_distribution_plot do
  auc_by_control_and_model = auc_by_control_and_model('occurences/auc/')
  max_aucs_for_uniprots = SequenceDataset.each_uniprot.map{|uniprot|
    SequenceDataset.each_for_uniprot(uniprot).map{|control|
      auc_by_control_and_model[control.name].map{|model, auc| auc }.max
    }.max
  }

  auc_distribution = (0.0..1.0).step(0.01).map{|t|
    cnt = max_aucs_for_uniprots.count{|x| x <= t }
    [t, cnt]
  }
  File.write 'auc_distribution.tsv', auc_distribution.map{|t, cnt| "#{t}\t#{cnt}" }.join("\n")
end
