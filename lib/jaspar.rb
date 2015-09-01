require 'uniprot_info'

module Jaspar
  MatrixInfo = Struct.new(:matrix_id, :collection, :base_id, :version, :name, :pcm_matrix, :uniprot_acs, :uniprot_ids, :species_ids, :species_names) do
    # Actually uniprot_ac can also be refseq / EMBL identifier
    def full_name; "#{base_id}.#{version} #{name}"; end
    def matrix_str(matrix_name: nil)
      matrix_wo_name_str = pcm_matrix.map{|row| row.join("\t") }.join("\n")
        ">#{matrix_name || full_name}\n#{matrix_wo_name_str}"
    end
    def to_s; full_name; end
    def inspect; to_s; end
    def fit_species?(species_name)
      species_names.include?(species_name)
    end
  end

  class Infos
    def initialize(position_counts_filename: 'MATRIX_DATA.txt', taxonomy_filename: 'TAX.txt',
                  matrix_species_filename: 'MATRIX_SPECIES.txt', matrix_proteins_filename: 'MATRIX_PROTEIN.txt',
                  matrix_name_infos_filename: 'MATRIX.txt', uniprot_infos:)
      @pcm_matrices_by_id = Infos.load_matrices(position_counts_filename)

      @species_id_by_species_name = Infos.load_species_id_by_name(taxonomy_filename)
      @species_name_by_species_id = Infos.load_species_name_by_id(taxonomy_filename)

      @species_ids_by_matrix_id = Infos.load_matrix_species(matrix_species_filename)
      # Actually uniprot_ac can also be refseq / EMBL identifier
      @matrix_uniprot_acs = Infos.load_matrix_uniprot_acs(matrix_proteins_filename)

      @matrix_name_info_by_id = MatrixNameInfo.each_in_file(matrix_name_infos_filename).map{|matrix_info|
        [matrix_info.matrix_id, matrix_info]
      }.to_h

      @uniprot_id_by_uniprot_ac = uniprot_infos.map{|info|
        [info.uniprot_ac, info.uniprot_id]
      }.to_h

      @uniprot_ids_by_refseq = group_list_by(uniprot_infos, :refseq_ids, :uniprot_id)
      @uniprot_ids_by_embl = group_list_by(uniprot_infos, :embl_ids, :uniprot_id)
    end

    def each_matrix
      return enum_for(:each_matrix)  unless block_given?
      @matrix_name_info_by_id.each_key do |matrix_id|
        matrix_name_info = @matrix_name_info_by_id[matrix_id]

        uniprot_acs = @matrix_uniprot_acs[matrix_id] || []
        uniprot_ids = uniprot_acs.flat_map{|uniprot_ac|
          # Actually uniprot_ac can also be refseq / embl
          [@uniprot_id_by_uniprot_ac[uniprot_ac], *@uniprot_ids_by_refseq[uniprot_ac], *@uniprot_ids_by_embl[uniprot_ac]]
        }.compact.uniq

        species_ids = @species_ids_by_matrix_id[matrix_id] || []
        species_names = species_ids.map{|species_id|  @species_name_by_species_id[species_id]  }.compact

        matrix_info = MatrixInfo.new(
          matrix_id, matrix_name_info.collection,
          matrix_name_info.base_id, matrix_name_info.version, matrix_name_info.name,
          @pcm_matrices_by_id[matrix_id],
          uniprot_acs, uniprot_ids,
          species_ids, species_names
        )

        yield matrix_info
      end
    end

    def group_list_by(list, groups_key, target_key)
      list.flat_map{|info|
        info.send(groups_key).map{|group|
          [group, info.send(target_key)]
        }
      }.group_by{|group, target|
        group
      }.map{|group, group_target_pairs|
        [group, group_target_pairs.map(&:last)]
      }.to_h
    end
    private :group_list_by

    def self.load_matrix_species(filename)
      File.readlines(filename).map{|line|
        matrix_id, species_ids = line.chomp.split("\t")
        [matrix_id, species_ids.split(',').map(&:strip)]
      }.to_h;
    end

    def self.load_species_id_by_name(filename)
      File.readlines(filename).map{|line|
        species_id, species_name = line.chomp.split("\t")
        [species_name, species_id]
      }.to_h
    end

    def self.load_species_name_by_id(filename)
      File.readlines(filename).map{|line|
        species_id, species_name = line.chomp.split("\t")
        [species_id, species_name]
      }.to_h
    end

    def self.load_matrix_uniprot_acs(filename)
      File.readlines(filename).map{|line|
        matrix_id, uniprot_acs = line.chomp.split("\t")
        [matrix_id, uniprot_acs.split(',').map(&:strip)]
      }.to_h;
    end


    # ID  row col val
    # 9229  A 1 0
    # 9229  A 2 3
    # ...
    def self.load_matrices(position_counts_filename)
      matrices = {}
      letter_indices = {'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3}
      File.readlines(position_counts_filename).drop(1).each{|line|
        matrix_id, letter, col, val = line.chomp.split("\t")
        matrices[matrix_id] ||= []
        pos = col.to_i - 1
        letter_index = letter_indices[letter]
        matrices[matrix_id][pos] ||= []
        matrices[matrix_id][pos][letter_index] = val.to_f
      }
      matrices
    end
  end

  # Additional data, representing matrix name and collection (CORE/CNE/PHYLOFACTS/SPLICE/POLII/FAM/PBM/PBM_HOMEO/PBM_HLH)
  MatrixNameInfo = Struct.new(:matrix_id, :collection, :base_id, :version, :name) do
    def self.each_in_file(filename, &block)
      File.readlines(filename).map{|line| self.from_string(line) }.each(&block)
    end

    def self.from_string(line)
      matrix_id, collection, base_id, version, name = line.chomp.split("\t")
      self.new(matrix_id, collection, base_id, version, name)
    end

    def to_s; "#{base_id}.#{version} #{name}"; end
    def inspect; to_s; end
  end
end
