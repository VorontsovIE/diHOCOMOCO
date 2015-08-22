def auc_by_control_and_model(folder)
  result = Hash.new{|h,k| h[k] = {} }
  Models.all_models.each{|model|
    model_auc_fn = File.join(folder, model.uniprot, "#{model.full_name}.txt")
    File.readlines(model_auc_fn).each{|line|
      control_name, auc = line.chomp.split("\t")
      result[control_name][model.full_name] = auc.to_f
    }
  }
  result
end

def auc_by_model_and_control(folder)
  result = Hash.new{|h,k| h[k] = {} }
  Models.all_models.each{|model|
    model_auc_fn = File.join(folder, model.uniprot, "#{model.full_name}.txt")
    File.readlines(model_auc_fn).each{|line|
      control_name, auc = line.chomp.split("\t")
      result[model.full_name][control_name] = auc.to_f
    }
  }
  result
end
