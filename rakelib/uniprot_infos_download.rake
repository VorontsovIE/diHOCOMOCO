desc 'Download Uniprot ID-AC mapping with lots of additional data'
task :uniprot_id_ac_mapping => 'uniprot_HomoSapiens_and_MusMusculus_lots_of_infos.tsv'

file 'uniprot_HomoSapiens_and_MusMusculus_lots_of_infos.tsv.gz' do
  query = 'organism:"Homo sapiens" OR organism:"Mus musculus"'
  columns = [
    'id', 'entry name', 'genes(PREFERRED)', 'genes(ALTERNATIVE)', 'protein names',
    'database(RefSeq)', 'database(EMBL)', 'organism-id',
    'database(HGNC)', 'database(MGI)', 'database(GeneID)', # GeneID is entrezgene
  ]

  options = {
    sort: 'score',
    desc: '',
    compress: 'yes',
    query: query,
    fil: '',
    format: 'tab',
    force: 'yes',
    columns: columns.join(','),
  }
  options_str = options.map{|k,v| "#{k}=#{v}" }.join('&')

  sh 'wget', "http://www.uniprot.org/uniprot/?#{options_str}", '-O', 'uniprot_HomoSapiens_and_MusMusculus_lots_of_infos.tsv.gz'
end

file 'uniprot_HomoSapiens_and_MusMusculus_lots_of_infos.tsv' => 'uniprot_HomoSapiens_and_MusMusculus_lots_of_infos.tsv.gz' do
  sh 'gzip', '--decompress', 'uniprot_HomoSapiens_and_MusMusculus_lots_of_infos.tsv.gz'
end
