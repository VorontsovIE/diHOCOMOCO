class TFClassification
  # terms by ids
  def initialize()
    @terms_by_id = {}
    @children_by_id = Hash.new{|h,k| h[k] = [] }
    @terms_by_name = Hash.new{|h,k| h[k] = [] }
    # @terms_by_id.each{|term_id, term|
    #   @children_by_id[term.parent_id] << term  if term.parent_id
    # }
  end

  def <<(term)
    raise "Duplicate id #{term.id}"  if @terms_by_id[term.id]
    @terms_by_id[term.id] = term
    @terms_by_name[term.name] << term
    @children_by_id[term.parent_id] << term  if term.parent_id
  end

  def self.from_file(filename)
    tf_ontology = self.new
    terms = File.readlines(filename)
      .map(&:chomp)
      .slice_before{|line|
        line.start_with?('[Term]')
      }.drop(1)
      .map{|enumerator|
        ClassificationTerm.from_line_array(tf_ontology, enumerator.to_a)
      }.each{|term|
        tf_ontology << term
      }
      tf_ontology
  end

  def term_by_name(name)
    @terms_by_name[name]
  end

  def term(term_id)
    @terms_by_id[term_id]
  end

  def children(term_id)
    @children_by_id[term_id]
  end

  def leaf?(term_id)
    raise "Term #{term_id} does not exist"  unless @terms_by_id[term_id]
    @children_by_id[term_id].empty?
  end

  def tf_groups(slice_deepness)
    @terms_by_id.each_value.select{|term|
      term.deepness >= slice_deepness && (!term.parent || term.parent.deepness < slice_deepness)
    }.map{|term|
      [term, term.subtree_nodes]
    }.to_h
  end

  ############

  ClassificationTerm = Struct.new(:ontology_tree, :id, :name, :subset, :definition, :parent_id, :uniprot_ACs, :other) do
    def self.from_line_array(ontology_tree, arr)
      id, name, subset, definition, parent_id = nil, nil, nil, nil, nil
      other = []
      uniprot_ACs = []

      arr.select{|line|
        line.match(/^\w+:/)
      }.each{|line|
        case line
        when /^id:/
          id = line[/^id: (?<data>.+)$/, :data]
        when /^name:/
          name = line[/^name: (?<data>.+)$/, :data]
        when /^subset:/
          subset = line[/^subset: (?<data>.+)$/, :data]
        when /^def:/
          definition = line[/^def: (?<data>.+)$/, :data]
        when /^is_a:/
          parent_id = line[/^is_a: (?<data>.+?) ! .+$/, :data]
        when /^xref: UNIPROT:/
          uniprot_ACs << line[/^xref: UNIPROT:(?<data>\w+)\b/, :data]
        else
          other << line
        end
      }

      self.new(ontology_tree, id, name, subset, definition, parent_id, uniprot_ACs, other)
    end

    def parent
      ontology_tree.term(parent_id)
    end

    def children
      ontology_tree.children(id)
    end

    def leaf?
      ontology_tree.leaf?(id)
    end

    # It can be different from number of ancestors
    def deepness
      id.split('.').size
    end

    def descendant_leafs
      leaf? ? [self] : children.flat_map(&:descendant_leafs)
    end

    def descendants
      children + children.flat_map(&:descendants)
    end

    def subtree_nodes
      [self] + children.flat_map(&:subtree_nodes)
    end

    def ancestors
      result = []
      term = self
      while term.parent
        term = term.parent
        result.unshift(term)
      end
      result
    end

    def to_s
      "#{name}{#{id}}"
    end

    def inspect; to_s; end
  end
end
