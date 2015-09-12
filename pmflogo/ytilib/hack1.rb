require 'rexml/formatters/pretty'

module REXML
  module Formatters
    # The Transitive formatter writes an XML document that parses to an
    # identical document as the source document.  This means that no extra
    # whitespace nodes are inserted, and whitespace within text nodes is
    # preserved.  Within these constraints, the document is pretty-printed,
    # with whitespace inserted into the metadata to introduce formatting.
    #
    # Note that this is only useful if the original XML is not already
    # formatted.  Since this formatter does not alter whitespace nodes, the
    # results of formatting already formatted XML will be odd.
    class Transitive < Default
      def initialize( indentation=2 )
        @indentation = indentation
        @level = 0
      end
      
      protected
      def write_element( node, output )
        output << "\n" << ' '*@level
        output << "<#{node.expanded_name}"

        node.attributes.each_attribute do |attr|
          output << " "
          attr.write( output )
        end unless node.attributes.empty?

        if node.children.empty?
          output << "/>"
        else
          output << ">"
          # If compact and all children are text, and if the formatted output
          # is less than the specified width, then try to print everything on
          # one line
          skip = false
          @level += @indentation
          
          only_text = true
          
          node.children.each { |child|
          	only_text = child.is_a?(REXML::Text) && only_text
            write( child, output )
          }
          @level -= @indentation
          output << "#{only_text ? "" : "\n" + ' '*@level}" << "</#{node.expanded_name}>"
        end
        
      end

      def write_text( node, output )
        output << node.to_s()
      end
    end
  end
  
  class Document
    def write( output=$stdout, indent=-1, trans=false, ie_hack=false )
      if xml_decl.encoding != "UTF-8" && !output.kind_of?(Output)
        output = Output.new( output, xml_decl.encoding )
      end
      formatter = if indent > -1
          if trans
            REXML::Formatters::Transitive.new( indent )
          else
            REXML::Formatters::Pretty.new( indent, ie_hack )
          end
        else
          REXML::Formatters::Default.new( ie_hack )
        end
      formatter.write( self, output )
    end
  end
end
