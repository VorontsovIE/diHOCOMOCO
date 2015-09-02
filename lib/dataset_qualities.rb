require 'rubyXL'

# similarity check: ["FAIL", "IMPOSSIBLE", "PASS", "-"]
# consistency check: ["IMPOSSIBLE", "NOT REQUIRED", "-", "FAIL", "PASS", "FAIl"]
# motif check: ["FAIL", "NOT REQUIRED", "WRONG", "PASS", "-", "CONSISTENT", "FAIl", "UNKNOWN"]
#
# chance to save: ["+", "NOT REQUIRED", "-"]
#
DatasetQuality = Struct.new(:dataset_name,
                            :di_selex_similarity, :di_hocomoco_similarity, :mono_selex_similarity, :mono_hocomoco_similarity,
                            :tf, :best_similarity,  :similarity_check, :consistency_check, :motif_check,  :chance_to_save) do

  def pass_quality_control?
    !([similarity_check, consistency_check, motif_check] & ['consistent', 'pass']).empty?
  end

  def self.nullable_float(x)
    x == 'NA' ? nil : x.to_f
  end

  def self.basename_wo_ext(filename)
    File.basename(filename, File.extname(filename))
  end

  def self.normalize_name(obj)
    obj.to_s.downcase
  end

  def model_names; [monoPWM_name, diPWM_name]; end
  def monoPWM_name; "#{tf}~CM~#{dataset_name}";  end
  def diPWM_name; "#{tf}~CD~#{dataset_name}";  end

  def models; [monoPWM_model, diPWM_model]; end
  def monoPWM_model; Model.new(monoPWM_name, :mono);  end
  def diPWM_model; Model.new(diPWM_name, :di);  end

  def self.from_row(row)
    chipseq_peak,  _diPWM, di_selex_similarity, di_hocomoco_similarity,  _monoPWM, mono_selex_similarity, mono_hocomoco_similarity,
    tf, best_similarity,  similarity_check, consistency_check, motif_check,  chance_to_save = row.size.times.map{|i| row[i] ? row[i].value : '' }
    self.new("#{tf}^#{chipseq_peak}", # dataset_name in table has no info about uniprot, so we extract its name from model name
            nullable_float(di_selex_similarity), nullable_float(di_hocomoco_similarity),
            nullable_float(mono_selex_similarity), nullable_float(mono_hocomoco_similarity),
            tf, nullable_float(best_similarity),
            normalize_name(similarity_check), normalize_name(consistency_check), normalize_name(motif_check),
            chance_to_save)
  end

  def self.each_in_xlsx(filename, worksheet: 0, &block)
    RubyXL::Parser.parse(filename)
      .worksheets[worksheet].drop(1)
      .map{|row|
        self.from_row(row)
      }.each(&block)
  end
end
