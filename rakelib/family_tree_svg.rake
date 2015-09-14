require 'models'
require 'best_models'
require 'auc_infos'
require 'motif_family_recognizer'

class TFClassification::ClassificationTerm
  def to_json(while_condition:, infos_proc:)
    return nil  unless while_condition.call(self)
    infos = infos_proc.call(self)
    result_params = infos.map{|k,v| "\"#{k}\":\"#{v}\""}

    children_nodes = children.sort.map{|child|
      child.to_json(while_condition: while_condition, infos_proc: infos_proc)
    }.compact
    result_params << ('"children":[' + children_nodes.join(",") + "]")  unless children_nodes.empty?
    "{" + result_params.map{|s| s }.join(",") + "}"
  end
end


def uniprots_in_subtree(node, uniprot_ids_by_ac, species)
  node.descendants
      .flat_map(&:uniprot_ACs)
      .flat_map{|uniprot_ac|
        uniprot_ids_by_ac[uniprot_ac]
      }.compact.uniq.select{|uniprot_id|
        uniprot_id.split('_').last == species
      }
  end

task :family_tree_svg do
  auc_infos_for_uniprot = AUCs.load_auc_infos_for_uniprot(min_weight_for_dataset: 0.65,
                                                          min_auc_for_model: 0.65);

  secondary_models = File.readlines('secondary_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) };
  banned_models = File.readlines('banned_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) };
  to_be_reversed = File.readlines('revcomp_models.txt').map(&:chomp).map{|name| Model.new_by_name(name) }.to_set;

  best_models = collect_best_models(auc_infos_for_uniprot,
                                    secondary_models: secondary_models,
                                    banned_models: banned_models);

  uniprot_infos = UniprotInfo.each_in_file('uniprot_HomoSapiens_and_MusMusculus_lots_of_infos.tsv');
  uniprot_ids_by_ac = uniprot_infos.group_by(&:uniprot_ac).map{|uniprot_ac, infos| [uniprot_ac, infos.map(&:uniprot_id)] }.to_h;


  ['HUMAN', 'MOUSE'].each do |species|
    ['mono', 'di'].each do |arity|

      best_models_subset = best_models.select{|model|
        model.species == species
      }.select{|model|
        model.arity_type == arity
      }
      uniprots_covered = best_models_subset.map(&:uniprot).uniq

      tf_ontology = TFClassification.from_file('TFOntologies/TFClass_human.obo')
      puts tf_ontology.root.to_json(
        while_condition: ->(node){
          uniprots_in_subtree = uniprots_in_subtree(node, uniprot_ids_by_ac, species)
          !uniprots_in_subtree.empty? && node.deepness <= 2
        },
        infos_proc:->(node){
          uniprots_in_subtree = uniprots_in_subtree(node, uniprot_ids_by_ac, species)

          uniprots_in_subtree_covered = uniprots_in_subtree.select{|uniprot_id|
            uniprots_covered.include?(uniprot_id)
          }

          {
            total_tfs: uniprots_in_subtree.size,
            size: uniprots_in_subtree.size,
            covered_tfs: uniprots_in_subtree_covered.size,
            name: node.name,
            family_id: node.id,
            url: "/search?species=#{species}&arity=#{arity}&family_id=#{node.id}"
          }
        }
      )
      exit
    end
  end
end

