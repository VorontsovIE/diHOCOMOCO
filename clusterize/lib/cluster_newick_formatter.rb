require_relative 'cluster_formatter'

class ClusterNewickFormatter < ClusterFormatter
  def create_newick_html(size={}, options = {})
    newick_inner_content = content()
    size_x = size[:x] || 8500
    size_y = size[:y] || 8500
    js_path = options[:js_path] || 'js'
    # css_path = options[:css_path] || 'css'
    
    <<-RESULT
    <html>
    <head>
      <script type="text/javascript" src="#{js_path}/raphael-min.js" ></script>
      <script type="text/javascript" src="#{js_path}/jsphylosvg-min.js"></script>
      
      <script type="text/javascript">
      window.onload = function(){
          var dataObject = { newick: "#{newick_inner_content};" };
          phylocanvas = new Smits.PhyloCanvas(
            dataObject,
            'svgCanvas', 
            #{size_x}, #{size_y},
            'circular'
          );
      };
      </script>
    </head>
    <body>
      <div id="svgCanvas"> </div>
    </body>
    </html>
    RESULT
  end

  def content_for_leaf_node(ind)
    compact_name( clusterer.names[ind] )
  end
  
  def content_for_inner_node(ind)
    ind1,ind2 = clusterer.children(ind)
    dist, dist1, dist2 = clusterer.send(branch_len_meth, ind), clusterer.send(branch_len_meth, ind1), clusterer.send(branch_len_meth, ind2)
    len_1 = (dist - dist1).round(3)
    len_2 = (dist - dist2).round(3)
    "(#{content_for_node(ind1)}:#{len_1},#{content_for_node(ind2)}:#{len_2})"
  end
  
  def content_for_node(ind)
    if clusterer.leaf?(ind)
      content_for_leaf_node(ind)
    else
      content_for_inner_node(ind)
    end
  end

  def content()
    content_for_node(clusterer.root_node)
  end
  
end