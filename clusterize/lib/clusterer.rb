require_relative 'cluster_newick_formatter'
require_relative 'cluster_xml_formatter'
require 'logger'
require 'yaml'

class Array
  # 2d-array dup
  def deep_dup
    dup.map!(&:dup)
  end
end

def load_matrix_from_file(matrix_filename)
  File.readlines(matrix_filename).map{|line| line.rstrip.split("\t").map(&:to_f) }
end

def load_matrix_from_file_with_names(matrix_filename)
  File.readlines(matrix_filename)[1..-1].map{|line| line.rstrip.split("\t")[1..-1].map(&:to_f) }
end

# Usage:
#  obj.eliminating_instance_variable(:@excess_data) { obj.dump }
class Object
  def eliminating_instance_variable(var)
    variable_backup = instance_variable_get(var)
    remove_instance_variable(var)
    yield
  ensure
    instance_variable_set(var, variable_backup)
  end
end

class Clusterer
  attr_accessor :names, :leafs_distance, :linkage_method
  attr_writer :logger, :max_distance

  def initialize(distance_matrix, linkage_method, names)
    raise ArgumentError, 'Negative distances in matrix'  if distance_matrix.any?{|line| line.any?{|el| el < 0.0 }}
    @size = distance_matrix.size
    raise 'Not squared matrix' unless distance_matrix.map(&:size).all?{|inner_sz| @size == inner_sz}
    @leafs_distance = distance_matrix
    @linkage_method = linkage_method
    @names = names
  end

  def indices_by_names(cluster_names)
    cluster_names.map{|name| names.index(name) }
  end

  def submatrix_by_indices(indices)
    leafs_distance.values_at(*indices).map{|row| row.values_at(*indices) }
  end

  def average_distances_in_cluster(cluster_names)
    return cluster_names.map{|name| [name, 0.0] }.to_h  if cluster_names.size <= 1
    inds = indices_by_names(cluster_names)
    mat = submatrix_by_indices(inds)
    cluster_names.zip(mat).map{|name, row_dists|
      [name, row_dists.inject(0.0, &:+) / (row_dists.size - 1)] # don't count distance of element to itself
    }.to_h
  end

  # find index (i,j) corresponding to a minimal element of distance submatrix(bounded on current_nodes).
  # return value ci,cj is always ordered: ci < cj
  # max_element is an initial value for min
  def index_of_minimal_element(current_nodes, current_dist, max_element)
    ci,cj = 0,1
    min = max_element
    current_nodes.each{|i|
      current_nodes.each{|j|
        next if i<=j
        if current_dist[i][j]  < min
          ci,cj = i,j
          min = current_dist[i][j]
        end
      }
    }
    return [ci,cj]
  end

  def update_distance_matrix(current_dist, ci, cj)
    sz = current_dist.size
    current_dist << sz.times.collect{|j| distance_to_cluster_union(current_dist, j, ci, cj) }
    sz.times{|i| current_dist[i] << distance_to_cluster_union(current_dist, i, ci, cj) } # duplicated calculations (not very time consuming)
    current_dist.last << 0.0
  end

  def make_linkage
    logger.info "Clustering started"
    tm = Time.now
    @linkage_tree = []
    current_nodes = num_items.times.to_a
    current_dist = leafs_distance #.deep_dup
    max_element = current_dist.map(&:max).max
    until current_nodes.size == 1
      ci, cj = index_of_minimal_element(current_nodes, current_dist, max_element)
      @linkage_tree << [[ci,cj], current_dist[ci][cj]]

      current_nodes << current_dist.size
      current_nodes.delete(cj)
      current_nodes.delete(ci)
      update_distance_matrix(current_dist, ci, cj)
      GC.start if current_dist.size % 100 == 0
      logger.info "#{current_dist.size - num_items} step -- #{Time.now - tm} sec"
    end
    logger.info "Clustering finished"
    self
  end

  def distance_to_cluster_union(distance_matrix, cx,ci,cj)
    send(linkage_method, distance_matrix, cx,ci,cj)
  end

  def single_linkage(distance_matrix, cx,ci,cj)
    [distance_matrix[cx][ci], distance_matrix[cx][cj]].min
  end
  def complete_linkage(distance_matrix,cx,ci,cj)
    [distance_matrix[cx][ci], distance_matrix[cx][cj]].max
  end
  def average_linkage(distance_matrix,cx,ci,cj)
    0.5 * (distance_matrix[cx][ci] + distance_matrix[cx][cj])
  end

  def num_items
    @size ||= leafs_distance.size
  end

  def root_node
    2 * num_items - 2
  end

  def leaf?(ind)
    ind < num_items
  end

  def children(ind) # return [ind1, ind2]
    @linkage_tree[ind - num_items].first
  end

  # glue clusters of subtrees when both subtrees aren't already divided into several clusters and if block returned true (or if block not given)
  def subtree_clusters(ind=root_node, &block)
    return [[ind]] if leaf?(ind)
    ind1,ind2 = children(ind)
    branch_1 = subtree_clusters(ind1,&block)
    branch_2 = subtree_clusters(ind2,&block)
    if branch_1.size == 1 && branch_2.size == 1 && (!block_given? || (block_given? && yield(ind)))
      [branch_1.first + branch_2.first]
    else
      branch_1 + branch_2
    end
  end

  def get_clusters_names(&block)
    subtree_clusters(&block).map{|clust| clust.map{|ind| names[ind]} }
  end

  def subtree_elements(ind)
    subtree_clusters(ind).flatten
  end

  def link_length(ind)
    leaf?(ind) ? 0.0 : @linkage_tree[ind - num_items].last
  end

  # max distance between all leaves in a subtree
  def subtree_max_distance(ind)
    #@max_distance[ind] ||= leaf?(ind) ? 0.0 : children(ind).map{|child| subtree_max_distance(child)}.max
    return max_distance[ind]  if max_distance[ind]
    nodes = subtree_elements(ind)
    max_distance = 0.0
    nodes.each{|i|
      nodes.each{|j|
        max_distance = (leafs_distance[i][j] > max_distance) ? leafs_distance[i][j] : max_distance
      }
    }
    @max_distance[ind] = max_distance
  end

  def inconsistences
    return @inconsistences  if @inconsistences
    result = []
    @linkage_tree.each do |((ind1,ind2),dst)|
      dst1 = link_length(ind1)
      dst2 = link_length(ind2)
      if dst1 == 0 && dst2 == 0
        result << -1
      else
        ave = (dst1 + dst2) / 2
        stddev = (dst2 - dst1).abs / Math.sqrt(2)
        result << (dst - ave) / stddev
      end
    end
    @inconsistences = result
  end

  def cutoff_criterium(criterium, cutoff)
    lambda{|ind| send(criterium, ind) <= cutoff}
  end

  def number_of_clusters_criterium(num_clusters)
    lambda{|ind| ind <= 2*num_items - 1 - num_clusters}
  end

  def dump_matrix(filename)
    File.open(filename,'w'){ |f|
      leafs_distance.each{|line| f.puts line.join("\t") }
    }
  end

  # dumps clustering tree without matrix (which can be dump separately if not yet saved by #dump_matrix) into YAML-file
  def dump(yaml_filename)
    eliminating_instance_variable(:@leafs_distance) {
      eliminating_instance_variable(:@logger) {
        File.open(yaml_filename,'w'){|f|  f.puts(self.to_yaml) }
      }
    }
  end

  def self.load(distance_matrix, yaml_filename)
    cluster = YAML.load_file(yaml_filename)
    cluster.leafs_distance = distance_matrix
    cluster
  end

  def max_distance
    @max_distance ||= []
  end

  def logger
    @logger ||= Logger.new($stderr)
  end

  def possible_cutoffs(criterium)
    (0..root_node).map{|ind| send(criterium,ind)}.uniq.sort
  end

  #  yields cutoff and array of clusters(names of motifs in cluster)
  #  at each possible cutoff
  def clusters_by_cutoff(criterium)
    return enum_for(:clusters_by_cutoff, criterium)  unless block_given?
    possible_cutoffs(criterium).each do |cutoff|
      clusters = get_clusters_names(&cutoff_criterium(criterium, cutoff))
      yield cutoff, clusters
    end
  end
end
