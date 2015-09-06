def print_html_table_for_grouped_models(auc_infos_for_uniprot, grouped_models, secondary_models:, stream: $stdout)
  stream.puts '<html><head><style>'
  stream.puts 'table, tr, td{ border: 1px solid black; border-collapse: collapse; }'
  stream.puts 'tr.mono, tr.mono td { border-top: 3px solid black; }'
  stream.puts 'tr.di, tr.di td{ border-bottom: 3px solid black; }'
  stream.puts 'tr.mono td:first-child{ border-bottom: 3px solid black; border-left: 3px solid;}'
  stream.puts 'tr td:last-child{ border-right: 3px solid; }'
  stream.puts '</style></head><body>'
  stream.puts '<table>'

  grouped_models.each do |uniprot, models|
    print_logos_for_models(auc_infos_for_uniprot[uniprot], models, uniprot, secondary_models: secondary_models, stream: stream)
  end

  stream.puts '</table>'
  stream.puts '</body></html>'
end

def print_logos_for_models(auc_infos, models, uniprot, secondary_models:, stream: $stdout)
  num_models = models.size
  models.each_with_index do |model, index|
    auc = auc_infos && auc_infos.weighted_model_aucs[model] && auc_infos.weighted_model_aucs[model].round(3)
    if index == 0
      stream.print '<tr class="mono">'
      stream.print "<td rowspan=#{num_models}>#{uniprot}</td>"
    else
      stream.print '<tr>'
    end

    quality = nil
    if secondary_models.include?(model)
      quality = 'S'
    end

    stream.print "<td>#{auc}</td>"
    stream.print "<td>#{quality}</td>"
    stream.print "<td>#{model.full_name}</td>"
    stream.print "<td><img src='#{model.path_to_logo_direct}'/></td>"
    stream.print "<td><img src='#{model.path_to_logo_revcomp}'/></td>"
    stream.puts '</tr>'
  end
end
