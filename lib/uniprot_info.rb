UniprotInfo = Struct.new(:uniprot_ac, :uniprot_id,
    :primary_gene_name, :synonym_gene_names, :protein_names,
    :refseq_ids, :embl_ids, :organism_id,
    :hgnc_ids, :mgi_ids, :entrezgene_ids) do
  def all_gene_names
    [primary_gene_name] + synonym_gene_names
  end

  def self.from_string(line)
    uniprot_ac, uniprot_id,  primary_gene_name, synonym_gene_names, \
      protein_names, refseq_ids, embl_ids, organism_id, \
      hgnc_ids, mgi_ids, entrezgene_ids = line.chomp.split("\t",11)

    refseq_ids = refseq_ids.split(';').map{|refseq_id| refseq_id.split('.').first }
    embl_ids = embl_ids.split(';')
    hgnc_ids = hgnc_ids.split(';')
    mgi_ids = mgi_ids.split(';')
    entrezgene_ids = entrezgene_ids.split(';')
    self.new(uniprot_ac, uniprot_id,
      primary_gene_name, synonym_gene_names.split,
      protein_names, refseq_ids, embl_ids, organism_id,
      hgnc_ids, mgi_ids, entrezgene_ids)
  end

  def self.each_in_file(filename, &block)
    File.readlines(filename).drop(1).map{|line| self.from_string(line) }.each(&block)
  end
end
