# This script should be run over final collection to check which di-motifs are oriented the wrong way
# Thrn final collection should be regenerated from scratch
task :which_di_to_reverse do
  File.open('curation/to_reverse_di.txt', 'w') do |fw|
    Dir.glob('final_collection/di/pwm/*.dpwm').each do |dipwm|
      uniprot = File.basename(dipwm).split('.').first
      main_mono = Dir.glob("final_collection/mono/pwm/#{uniprot}.H11MO.0.*.pwm").first
      cmd = 'java -cp /home/ilya/iogen_tools/macro-perfectos-ape/classes/artifacts/macro_perfectos_ape_jar/ape.jar ' + \
            ' ru.autosome.macroape.di.EvalSimilarity ' + \
            " #{main_mono}  #{dipwm} " +\
            ' --first-from-mono ' + \
            ' | grep -P "^OR\b" | cut -f2 '
      fw.puts File.basename(dipwm, '.dpwm')  if `#{cmd}`.strip == 'revcomp'
    end
  end
end
