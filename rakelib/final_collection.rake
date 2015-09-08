require 'models'
require 'best_models'
require 'auc_infos'
require 'quality_assessor'
require 'html_table_output'

desc 'Collect final collection'
task :make_final_collection do
  auc_infos_for_uniprot = AUCs.load_auc_infos_for_uniprot(min_weight_for_dataset: 0.65,
                                                          min_auc_for_model: 0.65)

  secondary_models = File.readlines('secondary_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }
  banned_models = File.readlines('banned_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }

  quality_assessor = QualityAssessor.new(auc_infos_for_uniprot, secondary_models)

  best_models = collect_best_models(auc_infos_for_uniprot,
                                    secondary_models: secondary_models,
                                    banned_models: banned_models)

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
  ['HUMAN', 'MOUSE'].each{|species|
    {'mono' => 'H10M', 'di' => 'H10D'}.each{|arity, bundle_name|
      mkdir_p "final_bundle/#{species}/#{arity}/pcm/"
      mkdir_p "final_bundle/#{species}/#{arity}/pwm/"
      mkdir_p "final_bundle/#{species}/#{arity}/logo/"

      best_models.select{|model|
        model.species == species
      }.select{|model|
        model.arity_type == arity
      }.group_by{|model|
        # All models from the same collection with the same original name refer to the same model
        # We will join these models into one
        [model.collection_short_name, model.model_name].join('~')
      }.each{|original_model_name, models|
        models = models.sort
        joint_uniprot = models.map(&:uniprot).join('+')

        qualities = models.map{|model| quality_assessor.calculate_quality(model) }
        if qualities.uniq.size != 1
          raise "There are different qualities for a model `#{model}` related to different TFs: `#{joint_uniprot}`" +
                '(the model was assigned to several UniprotIDs but copies got different qualities)'
        end
        quality = qualities.first

        model = models.first # All models are the same
        unless models.all?{|m| m.pcm.matrix == model.pcm.matrix && m.pwm.matrix == model.pwm.matrix }
          raise 'We expect all models with the same name to be the same'
        end
        model_name = "#{joint_uniprot}~#{bundle_name}~#{quality}"

        models_good_strand = models.map{|model| !to_be_reversed.include?(model) }
        if models_good_strand.uniq.size != 1
          raise "Model #{model_name} to be reversed and not be reversed simultaneousely"
        end
        good_strand = models_good_strand.first

        if good_strand
          pcm = model.pcm.named(model_name)
          pwm = model.pwm.named(model_name)
          logo_direct_path = model.path_to_logo_direct
          logo_revcomp_path = model.path_to_logo_revcomp
        else
          pcm = model.pcm.named(model_name).revcomp
          pwm = model.pwm.named(model_name).revcomp
          logo_direct_path = model.path_to_logo_revcomp
          logo_revcomp_path = model.path_to_logo_direct
        end

        folder = "final_bundle/#{model.species}/#{model.arity_type}"
        File.write File.join(folder, 'pcm', "#{model_name}.#{model.pcm_extension}"), pcm.to_s
        File.write File.join(folder, 'pwm', "#{model_name}.#{model.pwm_extension}"), pwm.to_s
        cp logo_direct_path, File.join(folder, 'logo', "#{model_name}_direct.png")
        cp logo_revcomp_path, File.join(folder, 'logo', "#{model_name}_revcomp.png")
      }
    }
  }
end
