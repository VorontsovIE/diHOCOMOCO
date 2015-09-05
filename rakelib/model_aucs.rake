require 'auc_infos'
require 'models'

task :model_logos_table do
  uniprots = FileList['occurences/auc/*'].pathmap('%n');
  auc_infos_for_uniprot = uniprots.map{|uniprot|
    [uniprot, AUCs.from_folder("occurences/auc/#{uniprot}/*.txt")]
  }.map{|uniprot, auc_infos|
    auc_infos_prev = nil
    while auc_infos != auc_infos_prev
      auc_infos_prev = auc_infos
      auc_infos = auc_infos.without_bad_datasets(0.65).without_bad_models(0.65)
    end
    [uniprot, auc_infos]
  }.reject{|uniprot,auc_infos|
    auc_infos.empty?
  }.to_h;

  puts '<html><head><style>'
  puts 'table, tr, td{ border: 1px solid black; border-collapse: collapse; }'
  puts 'tr.mono, tr.mono td { border-top: 3px solid black; }'
  puts 'tr.di, tr.di td{ border-bottom: 3px solid black; }'
  puts 'tr.mono td:first-child{ border-bottom: 3px solid black; border-left: 3px solid;}'
  puts 'tr td:last-child{ border-right: 3px solid; }'
  puts '</style></head><body>'
  puts '<table>'

  num_models = num_mono_models = num_di_models = 0
  auc_infos_for_uniprot.each{|uniprot, auc_infos|
    best_model_mono = auc_infos.best_model_among_collections(Models::CollectionsForFinalBundle & Models::MonoCollections)
    best_model_di = auc_infos.best_model_among_collections(Models::CollectionsForFinalBundle & Models::DiCollections)

    next  unless best_model_mono || best_model_di

    num_models += 1

    print '<tr class="mono">'
    print "<td rowspan=2>#{uniprot}</td>"
    if best_model_mono
      num_mono_models += 1
      best_auc_mono = auc_infos.weighted_model_aucs[best_model_mono]
      print "<td>#{best_auc_mono.round(3)}</td>"
      print "<td>#{best_model_mono.full_name}</td>"
      print "<td><img src='#{best_model_mono.path_to_logo}'/></td>"
    end
    print '</tr>'
    print '<tr class="di">'
    if best_model_di
      num_di_models += 1
      best_auc_di = auc_infos.weighted_model_aucs[best_model_di]
      print "<td>#{best_auc_di.round(3)}</td>"
      print "<td>#{best_model_di.full_name}</td>"
      print "<td><img src='#{best_model_di.path_to_logo}'/></td>"
    end
    puts "</tr>"
  }

  puts '</table>'
  puts '</body></html>'

  $stderr.puts({models: num_models, mono_models: num_mono_models, di_models: num_di_models})
end
