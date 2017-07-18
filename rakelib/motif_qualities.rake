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

task 'print_motif_qualities' do
  to_reverse = File.readlines('curation/revme_fin.txt').map(&:chomp)

  pcm_ext = {'mono' => 'pcm', 'di' => 'dpcm'}
  pwm_ext = {'mono' => 'pwm', 'di' => 'dpwm'}
  collection_name = {'mono' => 'H11MO', 'di' => 'H11DI'}

  ['mono', 'di'].each do |model_kind|
    FileUtils.mkdir_p "final_collection/#{model_kind}/pcm/"
    FileUtils.mkdir_p "final_collection/#{model_kind}/pwm/"
    FileUtils.mkdir_p "final_collection/#{model_kind}/logo/"
    File.open("final_collection/#{model_kind}.html", 'w') do |fw|
      fw.puts <<-EOS
<html><head></head><body><table>
      EOS

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
            if should_reverse
              if model_kind == 'mono'
                pcm = Bioinform::MotifModel::PCM.from_file("models/pcm/mono/all/#{model.split('~')[0]}/#{model}.pcm")
                pwm = Bioinform::MotifModel::PWM.from_file("models/pwm/mono/all/#{model.split('~')[0]}/#{model}.pwm")
              else
                pcm = ModelKind::Di.new.read_pcm("models/pcm/di/all/#{model.split('~')[0]}/#{model}.dpcm")
                pwm = ModelKind::Di.new.read_pwm("models/pwm/di/all/#{model.split('~')[0]}/#{model}.dpwm")
              end
              File.write("final_collection/#{model_kind}/pcm/#{final_name}.#{pcm_ext[model_kind]}", pcm.revcomp.to_s)
              File.write("final_collection/#{model_kind}/pwm/#{final_name}.#{pwm_ext[model_kind]}", pwm.revcomp.to_s)
            else
              FileUtils.cp "models/pcm/#{model_kind}/all/#{model.split('~')[0]}/#{model}.#{pcm_ext[model_kind]}", "final_collection/#{model_kind}/pcm/#{final_name}.#{pcm_ext[model_kind]}"
              FileUtils.cp "models/pwm/#{model_kind}/all/#{model.split('~')[0]}/#{model}.#{pwm_ext[model_kind]}", "final_collection/#{model_kind}/pwm/#{final_name}.#{pwm_ext[model_kind]}"
            end
            # puts [final_name, model].join("\t")
            fw.puts "<tr><td>#{final_name}</td><td>#{model}</td><td><img height='50' src='#{model_kind}/logo/#{final_name}_direct.png'></td></tr>"
          }
        }
      end
      fw.puts <<-EOS
</table></body></html>
      EOS


    end
  end
end
