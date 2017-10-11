DATASET_QUALITY_CUTOFF=$1
echo $'TF\tModel\tTPR at Pvalue=0.001\tTPR at Pvalue=0.0005\tTPR at Pvalue=0.0001';
for MOTIF in $(find final_bundle/hocomoco11/core/HUMAN/mono/pwm/ -xtype f | xargs -n1 basename -s .pwm | sort); do
	TF=$(ruby -e 'puts ARGV[0].split(".")[0]' -- $MOTIF)
	TF_WO_SPECIES=$( ruby -e 'puts ARGV[0].split("_")[0]' -- $TF )
	MOTIF_PATH="final_bundle/hocomoco11/core/HUMAN/mono/pwm/${MOTIF}.pwm"
	MODEL_LENGTH=$(( $(cat $MOTIF_PATH | wc -l) - 1 ))
	
	[ -z "$(find control/control/ -name "$TF.*.mfa")" ] && continue
	[ -z "$(find wauc_datasets/mono/ -name "${TF_WO_SPECIES}.*")" ] && continue
	echo -n $TF $'\t' $MOTIF $'\t';
	
	ls wauc_datasets/mono/${TF_WO_SPECIES}.* \
	  | ruby -e 'readlines.each{|fn| puts File.read(fn.chomp) }' \
	  | ruby -e "puts readlines.map{|l| nm, w = l.chomp.split(\"\\t\"); [nm,Float(w)] }.select{|nm,w| w >= $DATASET_QUALITY_CUTOFF }.map(&:first).uniq" \
	  | xargs -r -n1 -I {CNTRL} echo "cat control/control/{CNTRL}.mfa | ruby renameMultifastaSequences.rb all"  | bash \
	  | java -cp sarus.jar ru.autosome.SARUS - ${MOTIF_PATH} besthit \
	  | grep -v '>all' | cut -f1 \
	  | ruby -e "thrs = File.readlines('final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt').map{|l| l.chomp.split(\"\\t\") }.map{|m,*vals| [m, vals.map(&:to_f) ] }.detect{|m,vals| m == '$MOTIF'}.last; vals = readlines.map(&:to_f); puts thrs.map{|t| vals.count{|x| x >= t }.to_f / vals.size }.join(\"\\t\")"
done
