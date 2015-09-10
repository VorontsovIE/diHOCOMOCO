require 'quality_assessor'
require 'joint_model'

def cell_html(content, html_options = {})
  if html_options.empty?
    "<td>#{content}</td>"
  else
    options_str = html_options.map{|k,v| "#{k}=#{v}" }.join(' ') # doesn't add quotation marks
    "<td #{options_str}>#{content}</td>"
  end
end

def img_html(src)
  "<img src=\"#{src}\"/>"
end

def row_classes(index, num_total)
  classes = []
  classes << 'start' if index == 0
  classes << 'finish' if index + 1 == num_total
  classes << 'intermediate' if classes.empty?
  classes
end

def print_html_table(stream:, &block)
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
  yield
  stream.puts '</table>'
  stream.puts '</body></html>'
end

###############

def print_html_table_for_grouped_models(auc_infos_for_uniprot, grouped_models, quality_assessor, stream: $stdout)
  print_html_table(stream: stream) do
    grouped_models.each do |uniprot, models|
      print_logos_for_models(auc_infos_for_uniprot[uniprot], models, uniprot, quality_assessor, stream: stream)
    end
  end
end

def print_logos_for_models(auc_infos, models, uniprot, quality_assessor, stream: $stdout)
  num_models = models.size
  models.each_with_index do |model, index|
    auc = auc_infos.weighted_auc(model) && auc_infos.weighted_auc(model).round(3)
    quality = quality_assessor.calculate_quality(model)

    classes = row_classes(index, num_models)

    stream.print "<tr class=\"#{classes.join(' ')}\">"
    stream.print cell_html(uniprot, rowspan: num_models)  if index == 0
    stream.print cell_html(quality)
    stream.print cell_html(auc)
    stream.print cell_html(model.full_name)
    stream.print cell_html(img_html(model.path_to_logo_direct))
    stream.print cell_html(img_html(model.path_to_logo_revcomp))
    stream.puts '</tr>'
  end
end

################

def print_html_table_by_model_infos(model_infos, stream: $stdout)
  print_html_table(stream: stream) do
    model_infos.group_by(&:uniprot).each do |uniprot, model_infos_group|
      print_html_row_for_model_group(model_infos_group, stream: stream)
    end
  end
end

def print_html_row_for_model_group(model_infos, stream: $stdout)
  num_models = model_infos.size
  model_infos.each_with_index do |model_info, index|
    model_name = model_info.full_name
    auc = model_info.auc && model_info.auc.round(3)

    classes = row_classes(index, num_models)

    stream.print "<tr class=\"#{classes.join(' ')}\">"
    stream.print cell_html(model_info.uniprot, rowspan: num_models)  if index == 0
    stream.print cell_html(model_info.quality)
    stream.print cell_html(auc)
    stream.print cell_html(model_name)
    stream.print cell_html(img_html("logo/#{model_name}_direct.png"))
    stream.print cell_html(img_html("logo/#{model_name}_revcomp.png"))
    stream.print cell_html(model_info.origin_models.map(&:full_name).join(', '))
    stream.print cell_html(model_info.comments.join("<br/>"))
    stream.puts '</tr>'
  end
end

################

def print_csv_table(model_infos, stream: $stdout)
  model_infos.sort_by(&:full_name).each_with_index do |model_info|
    infos = [
      model_info.uniprot,
      model_info.quality,
      model_info.auc,
      model_info.full_name,
      model_info.origin_models.map(&:full_name).join(', '),
      model_info.comments.join(" "),
    ]
    stream.puts infos.join("\t")
  end
end
