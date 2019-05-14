require_relative 'cluster_formatter'

class ClusterXMLFormatter < ClusterFormatter
  attr_accessor :cutoff
  def initialize(clusterer, branch_len_meth, cutoff)
    super(clusterer, branch_len_meth)
    @cutoff = cutoff
  end

  def create_html_connected_to_xml(xml_file, size={}, options = {})
    size_x = size[:x] || 8500
    size_y = size[:y] || 8500
    js_path = options[:js_path] || 'js'
    css_path = options[:css_path] || 'css'

    <<-RESULT
    <html>
    <head>
      <script type="text/javascript" src="#{js_path}/yui-min.js"></script>
      <script type="text/javascript" src="#{js_path}/raphael-min.js" ></script>
      <script type="text/javascript" src="#{js_path}/jsphylosvg-min.js"></script>
      <script type="text/javascript" src="#{js_path}/unitip.js"></script>
      <link rel="stylesheet" type="text/css" href="#{css_path}/cssfonts-min.css">
      <link rel="stylesheet" type="text/css" href="#{css_path}/unitip.css">

      <script type="text/javascript">
      window.onload = function(){
        YUI().use('oop', 'json-stringify', 'io-base', 'event', 'event-delegate', function(Y){
          var uri = "#{xml_file}";
          function complete(id, o, args) {
            var data = o.responseXML; // Response data.
            var dataObject = {
                  xml: data,
                  fileSource: true
                };
            phylocanvas = new Smits.PhyloCanvas(
              dataObject,
              'svgCanvas',
              #{size_x}, #{size_y},
              'circular'
            );

            document.getElementById('svgCode').value = phylocanvas.getSvgSource();
            init(); //unitip
          };
          Y.on('io:complete', complete, Y);
          var request = Y.io(uri);
        });
      };
      </script>
    </head>
    <body>
      <div id="svgCanvas"> </div>
      <textarea id="svgCode"></textarea>
    </body>
    </html>
    RESULT
  end

  def clusters
    @clusters ||= clusterer.subtree_clusters(&clusterer.cutoff_criterium(branch_len_meth, cutoff))
  end

  def content
    xml_inner_content = content_for_node(clusterer.root_node)

    <<-RESULT
    <phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns="http://www.phyloxml.org">
    <phylogeny rooted="false">
    <render>
      <!-- <charts>
        <component type="binary" thickness="10" isInternal="true" />
      </charts> -->
      <styles>
        <first_cluster_group fill='#9F9' stroke='#FFF' />
        <second_cluster_group fill='#F99' stroke='#FFF' />
      </styles>
      <parameters>
        <circular>
            <innerCircleRadius>100</innerCircleRadius>
        </circular>
      </parameters>
    </render>
    <clade>
    <branch_length> 0.0 </branch_length>
    <clade>
    <branch_length> 0.0 </branch_length>
    #{ xml_inner_content }
    </clade>
    </clade>
    </phylogeny>
    </phyloxml>
    RESULT
  end

  def content_for_leaf_node(ind)
    stylename = clusters.find_index{|cluster| cluster.include? ind} % 2 == 0 ? 'first_cluster_group' : 'second_cluster_group'
    name = clusterer.names[ind]
    species = name.split('.')[0].split('_').last
    <<-RETSTRING
      <name bgStyle="#{stylename}">#{compact_name(name)}</name>
      <annotation>
      <desc>&lt;img src=\"http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/#{species}/mono/logo_small/#{name}_direct.png\" class=\"motif_desc\"/&gt;</desc>
      <uri>http://hocomoco11.autosome.ru/motif/#{name}</uri>
      </annotation>
      <chart>
        <component>#{stylename}</component>
      </chart>
    RETSTRING
  end

  def content_for_inner_node(ind)
    ind1, ind2 = clusterer.children(ind)
    dist, dist1, dist2 = clusterer.send(branch_len_meth, ind), clusterer.send(branch_len_meth, ind1), clusterer.send(branch_len_meth, ind2)
    len_1 = (dist - dist1).round(3)
    len_2 = (dist - dist2).round(3)
    <<-RETSTRING
      <clade>
      <branch_length>#{len_1}</branch_length>
      #{content_for_node(ind1)}
      </clade>
      <clade>
      <branch_length>#{len_2}</branch_length>
      #{content_for_node(ind2)}
      </clade>
    RETSTRING
  end

  # block: ind in linkage tree  -->  dist
  def content_for_node(ind)
    if clusterer.leaf?(ind)
      content_for_leaf_node(ind)
    else
      content_for_inner_node(ind)
    end
  end
end