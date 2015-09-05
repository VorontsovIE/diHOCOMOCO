$:.unshift './lib'
require 'cgi'
require 'rake'
require 'models'
require 'best_models'
require 'auc_info'
require 'auc_infos_filtering'
require 'quality_assessor'

filtering = AUCInfosFiltering.new(AUCInfo.load_all_infos)
filtering.remove_bad_datasets_and_models!(0.65)

quality_assessor = QualityAssessor.new(filtering)
# collection_perfomances = filtering.max_model_perfomances_collections_grouped
collection_perfomances = filtering.median_model_perfomances_collections_grouped


puts '<html><head></head><body>'
puts '<table>'
collection_perfomances.each do |uniprot, model_infos_by_collection|
  best_infos_di = best_model_infos(model_infos_by_collection, quality_assessor, Models::CollectionsForFinalBundle & Models::DiCollections)
  best_infos_mono = best_model_infos(model_infos_by_collection, quality_assessor, Models::CollectionsForFinalBundle & Models::MonoCollections)

  next  unless best_infos_mono[:model]
  model = Model.new(best_infos_mono[:model], :mono)

  print '<tr>'
  print "<td>#{uniprot}</td>"
  print "<td>#{best_infos_mono[:auc].round(2)}</td>"
  print "<td>#{best_infos_mono[:model]}</td>"
  print "<td><img src='#{CGI.escape(model.path_to_logo)}'/></td>"
  puts "</tr>"
end
puts '</table>'
puts '</body></html>'


# models = Models.all_models_by_uniprot('SP1_HUMAN').map(&:full_name)
# aucs = models.map{|model| [model, filtering.auc_infos.auc_by_model_and_dataset[model]] }.to_h
# datasets = aucs.values.flat_map(&:keys).uniq

# best_aucs_by_datasets = datasets.map{|ds|
#   aucs_for_ds = models.map{|model|
#     [model, aucs[model][ds]]
#   }.to_h
#   [ds, aucs_for_ds.max_by{|k,v| v}]
# }.to_h


# best_aucs_by_datasets.each{|ds, (model, best_auc)|
#   puts "#{ds}\t#{model}\t#{best_auc}"
# }

