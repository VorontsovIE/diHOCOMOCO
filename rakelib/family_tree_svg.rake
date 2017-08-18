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
  uniprot_infos = UniprotInfo.each_in_file('uniprot_HomoSapiens_and_MusMusculus_lots_of_infos.tsv');
  uniprot_ids_by_ac = uniprot_infos.group_by(&:uniprot_ac).map{|uniprot_ac, infos| [uniprot_ac, infos.map(&:uniprot_id)] }.to_h;

  ['HUMAN', 'MOUSE'].each do |species|
    tf_ontology = TFClassification.from_file("TFOntologies/TFClass_#{species.downcase}.obo")

    ['mono', 'di'].each do |arity|
      uniprots_covered = File.readlines("final_bundle/#{species}/#{arity}/final_collection.tsv").drop(1).map{|l|
        l.chomp.split("\t")[3]
      }

      json_str = tf_ontology.root.to_json(
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
      File.open("final_bundle/#{species}/#{arity}/family_tree_#{species}_#{arity}.json", 'w'){|fw|
        fw.puts json_str
      }
    end
  end
end

