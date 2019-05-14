#  ruby -r webrick -e "server = WEBrick::HTTPServer.new(Port: 2000, DocumentRoot: Dir.pwd); trap('INT'){server.shutdown}; server.start;" &

require_relative 'lib/clusterer'
require 'yaml'
require 'optparse'

def motif_in_transfac_format(matrix, name)
  "//\n" + \
  "NA\t#{name}\n" + \
  "P0\tA\tC\tG\tT\n" + \
  matrix.each_with_index.map{|pos,ind|
    "#{1 + ind}\t" + pos.join("\t") + "\n"
  }.join + 
  + "//"
end

def select_central_motif(motif_distances)
  motif_distances.sort_by{|motif, dist|
    uniprot, collection, rank, quality = motif.split('.') # hocomoco11 motif naming convention
    # two motifs always (and several motifs occasionally) have the same distance to others
    # in that case we want to select the most quality motif
    # (and if it isn't possible, we get alphabetically first one to make procedure stable)
    [dist, quality, Integer(rank), uniprot]
  }.first.first
end

options = {size: {x: 2500, y: 2500}, js_path: 'js', css_path: 'css'}
OptionParser.new{|cmd|
  cmd.banner = 'Usage: ruby clusterize.rb <distance_matrix.txt>'
  cmd.on('--newick [NEWICK_FILE]', 'generate newick file'){ |newick_file|
    options[:newick_file] = newick_file || 'cluster.newick'
  }
  cmd.on('--cluster-list [FILE]', 'generate list of clusters, collection of motifs in transfac format and cluster-TFs associations for MARA'){ |file|
    options[:clusters_list_file] = file || 'cluster_list.txt'
  }
  cmd.on('--xml [XML_FILE]', 'generate xml file'){ |xml_file|
    options[:xml_file] = xml_file || 'cluster.xml'
  }
  cmd.on('--yaml [YAML_FILE]', 'generate yaml file'){ |yaml_file|
    options[:yaml_file] = yaml_file || 'cluster.yaml'
  }
  cmd.on('--html [HTML_FILE]', 'generate html file'){ |html_file|
    options[:html_file] = html_file || 'cluster.html'
    options[:xml_file] ||= 'cluster.xml' # html is based on xml so we should generate it
    # Different order of --html and --xml options can give different results
  }
  cmd.on('--js-path JS_PATH', 'path for copying JS libs (jsPhyloSVG, Raphael, YUI, Unitip) necessary for html rendering'){ |js_path|
    options[:js_path] = File.absolute_path(js_path)
  }
  cmd.on('--css-path CSS_PATH', 'path for copying CSS stylesheet necessary for html rendering'){ |css_path|
    options[:css_path] = File.absolute_path(css_path)
  }
  cmd.on('--size SIZE', 'size of html-canvas in pixels in format: WIDTH,HEIGHT') {|size|
    options[:size][:x], options[:size][:y] = size.split(',', 2).map(&:to_i)
  }
  cmd.on('--mode MODE', 'use `distance` or `similarity` mode (similarity is in range 0..1 and is inverse of distance)') {|distance_mode|
    unless %w[distance similarity].include?(distance_mode)
      $stderr.puts 'Mode should be either distance or similarity'
      exit 1
    end
    options[:distance_mode] = distance_mode.to_sym
  }
  cmd.on('--separate-clusters THRESHOLD', 'separate clusters by defined threshold and colorize them in different colors') do |separate_threshold|
    options[:separate_threshold] = separate_threshold
  end
}.parse!

matrix_filename = ARGV.shift

raise 'distance matrix not specified'  unless matrix_filename
raise "#{matrix_filename} not exist"  unless File.exist?(matrix_filename)

# FileUtils.mkdir_p(File.join(output_folder,'html'))  unless Dir.exist? File.join(output_folder,'html')
# FileUtils.mkdir_p(File.dirname(cluster_dump_filename))  if cluster_dump_filename && ! Dir.exist?(File.dirname(cluster_dump_filename))


if options[:distance_mode] == :distance
  distance_matrix = load_matrix_from_file_with_names(matrix_filename)
else
  similarity_matrix = load_matrix_from_file_with_names(matrix_filename)
  distance_matrix = similarity_matrix.map{|l| l.map{|el| 1.0 - el } }
end

names = File.open(matrix_filename){|f|f.readline}.rstrip.split("\t")[1..-1]

clusterer = Clusterer.new(distance_matrix, :average_linkage, names)
clusterer.make_linkage
LINK_LENGTH_THRESHOLD = 0.95

if options[:clusters_list_file]
  clusters_with_centers = clusterer.get_clusters_names(&clusterer.cutoff_criterium(:link_length, LINK_LENGTH_THRESHOLD)).map{|cluster|
    distances_in_cluster = clusterer.average_distances_in_cluster(cluster)
    [select_central_motif(distances_in_cluster), cluster]
  }.sort_by(&:first).map{|cluster_center, motifs|
    [cluster_center, motifs.sort]
  }

  dn = File.dirname(options[:clusters_list_file])
  en = File.extname(options[:clusters_list_file])
  bn = File.basename(options[:clusters_list_file], en)
  File.open(File.join(dn, "#{bn}.txt"), 'w') do |fw|
    clusters_with_centers.each{|central_motif, motifs|
      species = central_motif.split('.')[0].split('_').last
      pcm_lines = File.readlines("final_bundle/hocomoco11/full/#{species}/mono/pcm/#{central_motif}.pcm").map(&:chomp)
      name = pcm_lines.first.gsub(/^>\s*/, '')
      mat = pcm_lines.drop(1).map{|l| l.split("\t").map(&:to_f) }
      fw.puts motif_in_transfac_format(mat, name)
    }
  end

  uniprot_infos = File.readlines('uniprot_HomoSapiens_and_MusMusculus_lots_of_infos.tsv').map{|l|
    uniprot_ac, uniprot_id, gene_name, alias_names, protein_name, refseq, embl, organism_id, hgnc, mgi, gene_id = l.chomp.split("\t", 11)
    infos = {
      uniprot_ac: uniprot_ac,
      uniprot_id: uniprot_id,
      gene_name: gene_name,
      alias_names: alias_names.split(' ').reject(&:empty?).join(','),
      protein_name: protein_name,
      refseq: refseq.split(';').reject(&:empty?).join(';'),
      embl: embl.split(';').reject(&:empty?).join(';'),
      organism_id: organism_id,
      hgnc: hgnc.split(';').reject(&:empty?).join(';'),
      mgi: mgi.split(';').reject(&:empty?).join(';'),
      gene_id: gene_id.split(';').reject(&:empty?).join(';'),
    }
    [uniprot_id, infos]
  }.to_h

  File.open(File.join(dn, "#{bn}_assoc.txt"), 'w') do |fw|
    clusters_with_centers.each{|central_motif, motifs|
      motif_infos = motifs.map{|motif|
        uniprot_id = motif.split('.').first
        infos = uniprot_infos[uniprot_id]

        [1, infos[:gene_id], infos[:gene_name], infos[:protein_name], infos[:uniprot_id], infos[:alias_names]].join(':')
      }
      fw.puts [central_motif, nil, *motif_infos].join("\t")
    }
  end
end

if options[:yaml_file]
  clusterer.dump(options[:yaml_file])
end

if options[:xml_file]
  xml_formatter = ClusterXMLFormatter.new(clusterer, :link_length, LINK_LENGTH_THRESHOLD)
  File.open(options[:xml_file],'w'){|f| f << xml_formatter.content() }
end

if options[:newick_file]
  newick_formatter = ClusterNewickFormatter.new(clusterer, :link_length)
  File.open(options[:newick_file],'w'){|f| f << newick_formatter.content()}
end

if options[:html_file]
  xml_formatter = ClusterXMLFormatter.new(clusterer, :link_length, LINK_LENGTH_THRESHOLD)
  newick_formatter = ClusterNewickFormatter.new(clusterer, :link_length)
  File.open(options[:html_file],'w'){|f| f << newick_formatter.create_newick_html(options[:size], js_path: options[:js_path]) }
  File.open(File.join(File.dirname(options[:html_file]), File.basename(options[:html_file], File.extname(options[:html_file])) + '_xml.html'),'w') do |f| 
    f << xml_formatter.create_html_connected_to_xml(options[:xml_file],
                                                    options[:size],
                                                    js_path: options[:js_path],
                                                    css_path: options[:css_path])
  end
  ClusterFormatter.copy_js_libs(options[:js_path])
  ClusterFormatter.copy_css_libs(options[:css_path])
end
