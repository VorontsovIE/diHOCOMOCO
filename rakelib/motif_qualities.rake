PCM_EXT = {'mono' => 'pcm', 'di' => 'dpcm'}
PWM_EXT = {'mono' => 'pwm', 'di' => 'dpwm'}

def peakset(dataset)
  dataset[/\.(PEAKS\d+)\./, 1]
end

def num_peaksets_scoring_better(auc_by_ds, auc_threshold)
  auc_by_ds.select{|ds, auc|
    auc >= auc_threshold
  }.map{|ds, auc|
    peakset(ds)
  }.uniq.count
end

def quality_mark(auc_by_ds, original_species)
  optimal_auc = 0.8
  minimal_auc = 0.65
  auc_by_origin_species_ds = auc_by_ds.select{|ds, auc|
    ds_species = ds.split('.').first.split('_').last
    ds_species == original_species
  }

  # At least one peak-set should pass threshold on the peak-set of original species
  # But if it does, we count number of peak-sets passed that threshold over all datasets (both HUMAN and MOUSE)
  any_origin_species_ds_pass_minimal = num_peaksets_scoring_better(auc_by_origin_species_ds, minimal_auc) >= 1
  num_datasets_pass_minimal_auc = any_origin_species_ds_pass_minimal ? num_peaksets_scoring_better(auc_by_ds, minimal_auc) : 0

  any_origin_species_ds_pass_optimal = num_peaksets_scoring_better(auc_by_origin_species_ds, optimal_auc) >= 1
  num_datasets_pass_optimal_auc = any_origin_species_ds_pass_optimal ? num_peaksets_scoring_better(auc_by_ds, optimal_auc) : 0


  auc_by_origin_species_ds.any?{|ds, auc| auc >= optimal_auc }

  if num_datasets_pass_optimal_auc >= 2
    'A'
  elsif num_datasets_pass_optimal_auc >= 1 && num_datasets_pass_minimal_auc >= 2
    'B'
  elsif num_datasets_pass_optimal_auc >= 1 || num_datasets_pass_minimal_auc >= 2
    'C'
  elsif num_datasets_pass_minimal_auc >= 1
    'D'
  else
    'E'
  end
end


def aucs_from_fn(filename)
  return []  unless File.exist?(filename)
  File.readlines(filename).map{|l|
    ds, auc, logauc = l.chomp.split("\t")
    [ds, Float(auc)]
  }
end

def put_motifs_to_final(model, final_name, model_kind, should_reverse)
  if !should_reverse
    FileUtils.cp "models/pcm/#{model_kind}/all/#{model.split('~')[0]}/#{model}.#{PCM_EXT[model_kind]}", "final_collection/#{model_kind}/pcm/#{final_name}.#{PCM_EXT[model_kind]}"
    FileUtils.cp "models/pwm/#{model_kind}/all/#{model.split('~')[0]}/#{model}.#{PWM_EXT[model_kind]}", "final_collection/#{model_kind}/pwm/#{final_name}.#{PWM_EXT[model_kind]}"
  else
    if model_kind == 'mono'
      pcm = Bioinform::MotifModel::PCM.from_file("models/pcm/mono/all/#{model.split('~')[0]}/#{model}.pcm")
      pwm = Bioinform::MotifModel::PWM.from_file("models/pwm/mono/all/#{model.split('~')[0]}/#{model}.pwm")
    else
      pcm = ModelKind::Di.new.read_pcm("models/pcm/di/all/#{model.split('~')[0]}/#{model}.dpcm")
      pwm = ModelKind::Di.new.read_pwm("models/pwm/di/all/#{model.split('~')[0]}/#{model}.dpwm")
    end
    File.write("final_collection/#{model_kind}/pcm/#{final_name}.#{PCM_EXT[model_kind]}", pcm.revcomp.to_s)
    File.write("final_collection/#{model_kind}/pwm/#{final_name}.#{PWM_EXT[model_kind]}", pwm.revcomp.to_s)
  end

end

task 'print_motif_qualities' do
  to_reverse = File.readlines('curation/revme_fin.txt').map(&:chomp)

  collection_name = {'mono' => 'H11MO', 'di' => 'H11DI'}

  ['mono', 'di'].each do |model_kind|
    FileUtils.mkdir_p "final_collection/#{model_kind}/pcm/"
    FileUtils.mkdir_p "final_collection/#{model_kind}/pwm/"
    FileUtils.mkdir_p "final_collection/#{model_kind}/logo/"
    table = []

    ['HUMAN', 'MOUSE'].each do |species|
      Dir.glob("wauc/#{model_kind}/#{species}/*.txt").sort.map{|slice_fn|

        model = File.readlines(slice_fn).map{|l|
          model, logauc, maxlogauc = l.chomp.split("\t")
          [model, Float(logauc)]
        }.max_by{|model, logauc|
          logauc
        }.first
        model_second_way = File.readlines(slice_fn).first.split("\t").first

        raise 'Inconsistent best model'  unless model == model_second_way # foolproof check

        auc_by_ds = aucs_from_fn("auc/#{model_kind}/HUMAN_datasets/#{model}.txt") + aucs_from_fn("auc/#{model_kind}/MOUSE_datasets/#{model}.txt")

        [File.basename(slice_fn, '.txt'), model, quality_mark(auc_by_ds, species)]
      }.group_by{|slice, model, quality|
        slice.split('.').first # semiuniprot
      }.each{|semiuniprot, slices|
        slices.sort_by{|slice, model, quality|
          slice_type = slice.split('.').last[0]
          ['M','S','T'].index(slice_type)
        }.each_with_index{|(slice, model, quality), ind|
          uniprot = "#{semiuniprot}_#{species}"
          final_name = "#{uniprot}.#{collection_name[model_kind]}.#{quality}#{ind}"
          should_reverse = to_reverse.include?(model.split('~').last)
          put_motifs_to_final(model, final_name, model_kind, should_reverse)
          table << ['novel', final_name, model, "<img src='#{model_kind}/logo/#{final_name}_direct.png'>",]
        }
      }
    end

    novel_uniprots = table.map{|final_name, model, img|
      final_name.split('.').first
    }.uniq

    hocomoco10_motifs_to_add = Dir.glob("hocomoco10/*/#{model_kind}/pcm/*").map{|pcm_fn|
      File.basename(pcm_fn, File.extname(pcm_fn))
    }.reject{|motif|
      uniprot = motif.split('.').first
      novel_uniprots.include?(uniprot)
    }

    hocomoco10_motifs_to_add.map{|motif|
      motif.split('.')
    }.group_by{|uniprot, old_collection, quality|
      uniprot
    }.each{|uniprot, models|
      species = uniprot.split('_').last
      uniprot_models = models.sort_by{|uniprot, old_collection, quality|
        ['A', 'B', 'C', 'D', 'S'].index(quality)
      }
      main_model_quality = uniprot_models.first.last
      uniprot_models.each_with_index{|(uniprot, old_collection, quality), ind|
        original_motif = "#{uniprot}.#{old_collection}.#{quality}"
        novel_quality = (quality != 'S') ? quality : main_model_quality

        novel_quality = 'D'  if uniprot.start_with?('GLI2_')

        hocomoco10_logo_fn = "hocomoco10/#{species}/#{model_kind}/logo_large/#{original_motif}_direct.png"
        novel_motif_name = "#{uniprot}.#{collection_name[model_kind]}.#{novel_quality}#{ind}"
        FileUtils.cp("hocomoco10/#{species}/#{model_kind}/logo_large/#{original_motif}_direct.png", "final_collection/#{model_kind}/logo/#{novel_motif_name}_direct.png")
        FileUtils.cp("hocomoco10/#{species}/#{model_kind}/logo_large/#{original_motif}_revcomp.png", "final_collection/#{model_kind}/logo/#{novel_motif_name}_revcomp.png")

        FileUtils.cp("hocomoco10/#{species}/#{model_kind}/pcm/#{original_motif}.#{PCM_EXT[model_kind]}", "final_collection/#{model_kind}/pcm/#{novel_motif_name}.#{PCM_EXT[model_kind]}")
        FileUtils.cp("hocomoco10/#{species}/#{model_kind}/pwm/#{original_motif}.#{PWM_EXT[model_kind]}", "final_collection/#{model_kind}/pwm/#{novel_motif_name}.#{PWM_EXT[model_kind]}")
        table << [
          'inherited',
          novel_motif_name,
          original_motif,
          "<img src='#{model_kind}/logo/#{novel_motif_name}_direct.png'>"
        ]
      }
    }

    motifs_to_ban = ['ERF', 'ETV2_HUMAN', 'MNT_HUMAN.H10MO.D', 'MUSC_HUMAN.H10MO.D', 'SMRC1', 'ZNF639', 'CLOCK_*.H10MO', 'PKNX2', 'YBX1', 'KAISO_MOUSE.H10MO.B']
    table.reject!{|inherit, final_name, model, img|
      uniprot = final_name.split('.').first
      motifs_to_ban.any?{|motif_pattern|
        pattern = /^(#{motif_pattern}\b|#{motif_pattern}_)/
        uniprot.match(pattern) || model.match(pattern)
      }
    }

    table.sort_by{|inherit, final_name, model, img|
      final_name
    }.chunk(&:itself).each_slice(300).map{|slice| slice.flat_map(&:last) }
    .each_with_index do |slice, slice_index|
      File.open("final_collection/#{model_kind}_slice_#{slice_index}.html", 'w'){|fw|
        fw.puts "<html><head><style>img{ height:50px; }\ntable,tr,td{ border:1px solid black; }\ntd:first-child{font-weight:bolder;}</style></head><body><table>"
        slice.each do |row|
          fw.puts("<tr>" + row.map{|cell| "<td>#{cell}</td>" }.join + "</tr>")
        end
        fw.puts "</table></body></html>"
      }
    end
  end
end
