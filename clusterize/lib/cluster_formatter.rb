require 'fileutils'

class ClusterFormatter
  attr_accessor :clusterer, :branch_len_meth
  def initialize(clusterer, branch_len_meth)
    @clusterer = clusterer
    @branch_len_meth = branch_len_meth
  end
  def compact_name(name)
    # mtch = /(?<rep>.+).+?\+\k<rep>/.match(name)
    # mtch ? name.gsub(/\+#{mtch[:rep]}/,'+') : name
    name.gsub(/\W/,'_')
  end
  def self.copy_js_libs(dest)
    FileUtils.cp_r(File.absolute_path('../../templates/js', __FILE__), dest)
  end
  def self.copy_css_libs(dest)
    FileUtils.cp_r(File.absolute_path('../../templates/css', __FILE__), dest)
  end
end
