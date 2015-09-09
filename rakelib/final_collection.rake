require 'models'
require 'best_models'
require 'joint_model'
require 'auc_infos'
require 'quality_assessor'
require 'html_table_output'

desc 'Collect final collection'
task :make_final_collection do
  auc_infos_for_uniprot = AUCs.load_auc_infos_for_uniprot(min_weight_for_dataset: 0.65,
                                                          min_auc_for_model: 0.65)

  secondary_models = File.readlines('secondary_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }
  banned_models = File.readlines('banned_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }


  best_models = collect_best_models(auc_infos_for_uniprot,
                                    secondary_models: secondary_models,
                                    banned_models: banned_models)

  quality_assessor = QualityAssessor.new(auc_infos_for_uniprot, best_models: best_models, secondary_models: secondary_models)

  File.open('final_collection.html', 'w') do |fw|
    print_html_table_for_grouped_models(
      auc_infos_for_uniprot,
      best_models.group_by(&:uniprot),
      quality_assessor,
      stream: fw
    )
  end


  to_be_reversed = File.readlines('revcomp_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }.to_set

  rm_rf 'final_bundle'
  ['HUMAN', 'MOUSE'].each do |species|
    ['mono', 'di'].each do |arity|
      models = best_models.select{|model|
        model.species == species
      }.select{|model|
        model.arity_type == arity
      }

      # combine same models for several Tfs into joint models
      model_infos = JointModel.grouped_models_from_scratch(models, auc_infos_for_uniprot, quality_assessor, to_be_reversed)

      folder = "final_bundle/#{species}/#{arity}"
      mkdir_p File.join(folder, 'pcm')
      mkdir_p File.join(folder, 'pwm')
      mkdir_p File.join(folder, 'logo')

      File.open(File.join(folder, "final_collection.html"), 'w') do |fw|
        print_html_table_by_model_infos(model_infos, stream: fw)
      end

      model_infos.each do |model_info|
        model_info.save_model_pack_into_folder!(folder)
      end
    end
  end
end
