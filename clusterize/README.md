TODO: names_filename shouldn't be obligatory in case when distance matrix has filenames!

    ruby clusterize.rb  distance_matrix/distance_matrix_with_names.txt  distance_matrix/clustering_results  cluster.yaml


    Options:
    ruby clusterize.rb [--with-names] [--log log_file.log] <matrix-txt-file> <motif-names-yaml-file>
      <output-folder> [cluster-yaml-dump]

cluster yaml dump used in a pair of ways:
 - when specified dump exists - it's loaded (in such case time-consuming clusterization stage's eliminated)
 - when dump not exists - built clusterer'll be dumped to specified file (so that next time clusterization run, it could immediately get clusterization tree)
 
