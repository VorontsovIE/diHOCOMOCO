require 'quality_assessor'

def print_html_table_for_grouped_models(auc_infos_for_uniprot, grouped_models, quality_assessor, stream: $stdout)
  stream.puts '<html><head><style>'
  stream.puts 'table, tr, td{ border: 1px solid black; border-collapse: collapse; padding: 2px 4px; }'
  stream.puts 'tr.start td:first-child{ border-top: 3px solid black; border-bottom: 3px solid black; border-left: 3px solid black;}'
  stream.puts 'tr.start, tr.start td { border-top: 3px solid black; }'
  stream.puts 'tr.start td { padding-top: 5px; }'
  stream.puts 'tr.finish td { padding-bottom: 5px; }'
  stream.puts 'tr td:last-child{ border-right: 3px solid; }'
  stream.puts 'tr.finish, tr.finish td{ border-bottom: 3px solid black;}'
  stream.puts '</style></head><body>'
  stream.puts '<table>'

  grouped_models.each do |uniprot, models|
    print_logos_for_models(auc_infos_for_uniprot[uniprot], models, uniprot, quality_assessor, stream: stream)
  end

  stream.puts '</table>'
  stream.puts '</body></html>'
end

def print_logos_for_models(auc_infos, models, uniprot, quality_assessor, stream: $stdout)
  num_models = models.size
  models.each_with_index do |model, index|
    auc = auc_infos.weighted_auc(model) && auc_infos.weighted_auc(model).round(3)
    quality = quality_assessor.calculate_quality(model)

    if index == 0
      if num_models == 1
        stream.print '<tr class="start finish">'
      else
        stream.print '<tr class="start">'
      end
      stream.print "<td rowspan=#{num_models}>#{uniprot}</td>"
    elsif index + 1 == num_models
      stream.print '<tr class="finish">'
    else
      stream.print '<tr class="intermediate">'
    end

    stream.print "<td>#{quality}</td>"
    stream.print "<td>#{auc}</td>"
    stream.print "<td>#{model.full_name}</td>"
    stream.print "<td><img src=\"#{model.path_to_logo_direct}\"/></td>"
    stream.print "<td><img src=\"#{model.path_to_logo_revcomp}\"/></td>"
    stream.puts '</tr>'
  end
end


def print_html_table_by_model_infos(model_infos, stream: $stdout)
  stream.puts '<html><head><style>'
  stream.puts 'table, tr, td{ border: 1px solid black; border-collapse: collapse; padding: 2px 4px; }'
  stream.puts 'tr.start td:first-child{ border-top: 3px solid black; border-bottom: 3px solid black; border-left: 3px solid black;}'
  stream.puts 'tr.start, tr.start td { border-top: 3px solid black; }'
  stream.puts 'tr.start td { padding-top: 5px; }'
  stream.puts 'tr.finish td { padding-bottom: 5px; }'
  stream.puts 'tr td:last-child{ border-right: 3px solid; }'
  stream.puts 'tr.finish, tr.finish td{ border-bottom: 3px solid black;}'
  stream.puts '</style></head><body>'
  stream.puts '<table>'

  model_infos.group_by{|info|
    info[:uniprot]
  }.each do |uniprot, models|
    print_html_row_for_model_group(models, stream: stream)
  end

  stream.puts '</table>'
  stream.puts '</body></html>'
end

def print_html_row_for_model_group(model_infos, stream: $stdout)
  num_models = model_infos.size
  model_infos.each_with_index do |model_info, index|
    model_name = model_info[:model_name]
    uniprot = model_info[:uniprot]
    quality = model_info[:quality]
    auc = model_info[:auc] && model_info[:auc].round(3)
    comment = model_info[:comment]

    if index == 0
      if num_models == 1
        stream.print '<tr class="start finish">'
      else
        stream.print '<tr class="start">'
      end
      stream.print "<td rowspan=#{num_models}>#{uniprot}</td>"
    elsif index + 1 == num_models
      stream.print '<tr class="finish">'
    else
      stream.print '<tr class="intermediate">'
    end

    stream.print "<td>#{quality}</td>"
    stream.print "<td>#{auc}</td>"
    stream.print "<td>#{model_name}</td>"
    stream.print "<td><img src=\"logo/#{model_name}_direct.png\"/></td>"
    stream.print "<td><img src=\"logo/#{model_name}_revcomp.png\"/></td>"
    stream.print "<td>#{model_info[:origin_models].map(&:full_name).join(' ')}</td>"
    stream.print "<td>#{comment}</td>"
    stream.puts '</tr>'
  end
end
