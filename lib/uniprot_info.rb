UniprotInfo = Struct.new(:uniprot_ac, :uniprot_id,  :primary_gene_name, :synonym_gene_names, :protein_names, :organism_id) do
  def all_gene_names
    primary_gene_name + synonym_gene_names
  end

  def self.from_string(line)
    uniprot_ac, uniprot_id,  primary_gene_name, synonym_gene_names, protein_names, organism_id = line.chomp.split("\t")
    self.new(uniprot_ac, uniprot_id,  primary_gene_name, synonym_gene_names.split, protein_names, organism_id)
  end

  def self.each_in_file(filename, &block)
    File.readlines(filename).drop(1).map{|line| self.from_string(line) }.each(&block)
  end
end
