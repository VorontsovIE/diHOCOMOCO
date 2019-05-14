# ruby clusterize.rb distance_matrix/distance_matrix_with_names.txt  distance_matrix/clustering_results  cluster.yaml
#
# ruby clusterize.rb [--with-names] [--log log_file.log] <matrix-txt-file> <motif-names-yaml-file> <output-folder> [cluster-yaml-dump]
# cluster yaml dump used in a pair of ways:
#  - when specified dump exists - it's loaded (in such case time-consuming clusterization stage's eliminated)
#  - when dump not exists - built clusterer'll be dumped to specified file (so that next time clusterization run, it could immediately get clusterization tree)

require_relative 'lib/clusterer'
require 'yaml'
require 'fileutils'
require 'optparse'
require '/home/ilya/iogen/fantom_known_recalc/naming_conversions'

class Clusterer
  def num_annotated_by_possible_cutoffs(criterium)
    annotated = ->(cluster){ cluster.any?{|motif| $naming_conversions.known?(motif) } && cluster.any?{|motif| $naming_conversions.denovo?(motif) } }
    result = []
    clusters_by_cutoff(criterium) do |cutoff, clusters|
      annotated_clusters = clusters.select(&annotated)
      result << [cutoff, clusters.size, annotated_clusters.size]
    end
    result
  end

  def statistics_by_swissregulon_possible_cutoffs(criterium)
    result = []
    clusters_by_cutoff(criterium) do |cutoff, clusters|
      num_swissregulons = clusters.map{|cluster| cluster.count{|motif| $naming_conversions.swissregulon?(motif) } }
      not_null_swissregulons = num_swissregulons.select{|cnt| cnt > 0}
      several_swissregulons = num_swissregulons.select{|cnt| cnt > 1}
      result << [cutoff,
                 clusters.size,
                 not_null_swissregulons.size,
                 num_swissregulons.count{|cnt| cnt > 1},
                 not_null_swissregulons.inject(&:+).to_f / not_null_swissregulons.size,
                 several_swissregulons.inject(&:+).to_f / several_swissregulons.size
                ]
    end
    result
  end

  def swissregulon_cutoffs(criterium)
    cutoff_before_1, _clusters_before_1 = clusters_by_cutoff(criterium).slice_before do |cutoff, clusters|
      num_swissregulons = clusters.map{|cluster| cluster.count{|motif| $naming_conversions.swissregulon?(motif)} }
      num_swissregulons.any?{|size_sw_inside| size_sw_inside > 1}
    end.first.last
    cutoff_before_2, _clusters_before_2 = clusters_by_cutoff(criterium).slice_before do |cutoff, clusters|
      num_swissregulons = clusters.map{|cluster| cluster.count{|motif| $naming_conversions.swissregulon?(motif)} }
      num_swissregulons.count{|size_sw_inside| size_sw_inside > 1} >= 2
    end.first.last
    cutoff_before_3, _clusters_before_3 = clusters_by_cutoff(criterium).slice_before do |cutoff, clusters|
      num_swissregulons = clusters.map{|cluster| cluster.count{|motif| $naming_conversions.swissregulon?(motif) } }
      num_swissregulons.any?{|size_sw_inside| size_sw_inside > 2}
    end.first.last
    [cutoff_before_1, cutoff_before_2, cutoff_before_3]
  end
end

# Returns sorted list of cutoffs with number of clusters at this cutoff and number of annotated clusters
# Format: [[cutoff_1, num_clusters_1, num_annotated_clusters_1], ...]
# Selection condition should be specified with a block (because there are different known-denovo recognition methods)
# Block gets cluster (names of motifs in cluster) and should return true iff cluster is annotated
def calculate_statistics_by_possible_cutoffs(clusterer, criterium, &block)
  result = []
  clusterer.possible_cutoffs(criterium).each{|cutoff|
    clusters = clusterer.get_clusters_names(&clusterer.cutoff_criterium(criterium, cutoff))
    annotated_clusters = clusters.select(&block)
    result << [cutoff, clusters.size, annotated_clusters.size]
  }
  File.open("#{criterium}_cutoffs.txt",'w'){|f| f.puts result.map{|cutoff_stat| cutoff_stat.join("\t")}}
  result
end

# Looks for a cutoff that yields maximal number of annotated clusters (consisting of both known and denovo motif).
# Returns range from least and highest possible cutoff (prefer more broad or compact clusters)
def best_cutoff(clusterer, criterium, &block)
  statistics = calculate_statistics_by_possible_cutoffs(clusterer, criterium, &block)
  max_num_annotated_clusters = statistics.map{|cutoff, num_clusters, num_annotated_clusters| num_annotated_clusters }.max
  best_cutoffs_infos = statistics.select{|cutoff, num_clusters, num_annotated_clusters| num_annotated_clusters == max_num_annotated_clusters }
  best_cutoffs_infos.first.first .. best_cutoffs_infos.last.first
end


$naming_conversions = NamingConversions.new
options = { }
OptionParser.new{|cmd|
  cmd.banner = 'Usage: ruby clusterize.rb <distance_matrix/distance_matrix_with_names.txt>  <distance_matrix/clustering_results>  <distance_matrix/cluster.yaml>'
  cmd.on('-l', '--log LOG_FILE', 'log-file of clusterization process (by default stderr used)'){ |log_file|
    options[:log_file] = log_file
  }
  cmd.on('-n','--names FILENAME','load names from separate file (not from distance matrix first row/column)'){ |names_filename|
    options[:names_file] = names_filename  # distance_matrix/motifs_order.yaml
    raise "names file #{names_filename} not exist"  unless File.exist?(names_filename)
  }
}.parse!

matrix_filename = ARGV.shift        # distance_matrix/distance_macroape.txt
output_folder = ARGV.shift          # distance_matrix/clustering_results
cluster_dump_filename = ARGV.shift  # distance_matrix/cluster.yaml

raise 'matrix filename not specified'  unless matrix_filename
raise "matrix file #{matrix_filename} not exist"  unless File.exist?(matrix_filename)
raise 'output folder not specified'  unless output_folder

FileUtils.mkdir_p(output_folder)  unless Dir.exist? output_folder
FileUtils.mkdir_p(File.join(output_folder,'html'))  unless Dir.exist? File.join(output_folder,'html')
FileUtils.mkdir_p(File.dirname(cluster_dump_filename))  if cluster_dump_filename && ! Dir.exist?(File.dirname(cluster_dump_filename))


if options[:names_file]
  distance_matrix = load_matrix_from_file(matrix_filename)
  names = YAML.load_file(options[:names_file])
else
  distance_matrix = load_matrix_from_file_with_names(matrix_filename)
  names = File.open(matrix_filename){|f|f.readline}.rstrip.split("\t")[1..-1]
end

if cluster_dump_filename
  if File.exist?(cluster_dump_filename)
    clusterer = Clusterer.load(distance_matrix, cluster_dump_filename)
  else
    clusterer = Clusterer.new(distance_matrix, :average_linkage, names)
    clusterer.logger = Logger.new(options[:log_file])  if options[:log_file]
    clusterer.make_linkage
    clusterer.dump(cluster_dump_filename)
  end
else
  clusterer = Clusterer.new(distance_matrix, :average_linkage, names)
  clusterer.make_linkage
end

newick_formatter = ClusterNewickFormatter.new(clusterer, :link_length)
xml_formatter = ClusterXMLFormatter.new(clusterer, :link_length, 0.1)
File.open("#{output_folder}/html/macroape_linklength.html",'w'){|f| f << newick_formatter.create_newick_html()}
File.open("#{output_folder}/html/macroape_linklength.newick",'w'){|f| f << newick_formatter.content()}
File.open("#{output_folder}/html/macroape_linklength.xml",'w'){|f| f << xml_formatter.content()}
#File.open("#{output_folder}/html/macroape_linklength_w_xml.html",'w'){|f| f << xml_formatter.create_html_connected_to_xml("#{output_folder}/macroape_linklength.xml") }
ClusterFormatter.copy_js_libs("#{output_folder}/html/")


# # threshold found with criterium that most(just before first and second gluing-events ) swissregulons hadn't glued in clusters:
# # link_length 0.744384589  497     182     8       1.043956044     2 - just before a point where more than 2 swissregulons becomes glued in one cluster
# # link_length 0.667869685  616     189     1       1.005291005     2 - just before point where second cluster with two swissregulons appears.
# # distance_macroape_cutoff_linklength_grid = [0.744384589, 0.667869685]
# # This results are obtained by manual processing of results of following code:
# #
# [:link_length, :subtree_max_distance].each do |criterium|
#   clusterer.statistics_by_swissregulon_possible_cutoffs(criterium).tap{|results|
#     File.write "statistics_by_swissregulon(#{criterium}).txt", results.map{|res| res.join("\t") }.join("\n")
#   }
# end
# swissregulon_cutoffs = clusterer.swissregulon_cutoffs(:link_length)
# p swissregulon_cutoffs


#
# Cutoffs choosen by criterium that maximal number of clusters (on whole set known+lexicon+denovo) are annotated
# if several cutoffs exist, one chooses cutoff corresponding to most compactified clusters:
# subtree_max_distance  -->  distance_macroape_cutoff_grid= [0.9746]
# link_length  -->  distance_macroape_cutoff_grid = [0.9586]
#
# Recalc gave us:
# subtree_max_distance --> 0.9995135346164128..0.9995430451218578
# link_length --> 0.9452718093750845..0.9586281289091381
#

annotated = ->(cluster){ cluster.any?{|motif| motif =~ /^KNOWN/} && cluster.any?{|motif| motif =~ /^DENOVO/ } }

distance_macroape_cutoff_grid = { subtree_max_distance: [], link_length: [] }

distance_macroape_cutoff_grid = { subtree_max_distance: [0.9995135346164128, 0.9995430451218578],
        link_length: [0.9452718093750845,0.9586281289091381] + swissregulon_cutoffs }

[:subtree_max_distance, :link_length].each do |criterium|
  clusters_macroape = {}
#  cutoffs = best_cutoff(clusterer, criterium, &annotated)
#  puts "criterium: #{criterium} -- #{cutoffs}"
#   cutoffs for most-compact clusters among clusters with maximal number of annotated ones
#  distance_macroape_cutoff_grid[criterium] << cutoffs.end
#   cutoffs for least-compact clusters among clusters with maximal number of annotated ones
#  distance_macroape_cutoff_grid[criterium] << cutoffs.begin

  distance_macroape_cutoff_grid[criterium].each do |cutoff|
    cutoff_criterium = clusterer.cutoff_criterium(criterium, cutoff)
    clusters_macroape[cutoff] = clusterer.get_clusters_names(&cutoff_criterium)
  end

  clusters_macroape.each do |cutoff, clusters|
    filename = "macroape_#{criterium}_cluster_names (#{cutoff.round(4)} - #{clusters.size}).txt"
    File.open("#{output_folder}/#{filename}",'w'){|f|
      clusters.each{|clust| f.puts clust.join("\t") }
    }
  end
end
