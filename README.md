# diHOCOMOCO
This code was used to prepare [https://github.com/VorontsovIE/diHOCOMOCO/releases/tag/beforeHocomoco11](Hocomoco-v10) and [https://github.com/VorontsovIE/diHOCOMOCO/releases/tag/hocomoco11](Hocomoco-v11) motif collection.
Its purpose is to assess and choose the best motif across multiple motifs. All the stages are described in Rakefile but it is not recommended to run it, it won't do what you can expect.

General workflow is the following:
* put PCM motifs, alignments and ChIP-seq control FASTA into certain folders. Please, prepare data on this stage manually.
Existing scripts heavily rely on data layout of previous steps.
* prepare PWMs from PCMs
* calculate average local dinucleotide background frequencies in ChIP-seq datasets of the corresponding TF
* estimate score-to-pvalue conversion for each model using calculated local background
* Calculate ROC AUC and logROC AUC values for each TF motif variant and each TF dataset. Cross-species experiments are also performed.
* Aggregate logROC AUC values into weighted logROC AUC (with weights corresponding to estimated dataset qualities). 
Motif packs (slices) to be aggregated are selected by manual curration (see motif list slices in `curation` folder).
Each slice corresponds to similar patterns of the same TF. 
Most TFs have a single slice, but sometimes there exist different binding modes (e.g. mono-box and double-box binding).
Datasets of secondary species are also taken into account with lower weight (in a sort of pseudocount addition procedure).
* Choose the best motif for each motif slice - motif with maximal weighted logROC AUC value. Collect all the data (pcm, pwm, thresholds, alignments) for these motifs.

Folder layout is the following:
* Models should be put into `models/pcm/(mono or di)/all/{TF}/` folder.
* Controls should be put into `control/control/` folder.
* Alignments go to `models/words/mono/` folders (in subfolders `hocomoco_legacy` and `chipseq`).

Model naming is the following:
* `{TF name}~{collection}~{motif name}.pcm` for mononucleoide motifs or `{TF name}~{collection}~{motif name}.dpcm` for dinucleotide motifs.
E.g. `KLF4_HUMAN~CM~KLF4_HUMAN.PEAKS034751.pics.S.pcm`
* Collection is one of the following abbreviations: 
`CM` (novel ChIP-seq beasd mononucleotide motifs), 
`CD` (novel ChIP-seq based dinucleotide motifs), 
`HL` (Hocomoco10 legacy motifs), 
`DIHL` (dinucleotide Hocomoco10 legacy motifs) and some others.

ChIP-seq control datasets are multi-FASTA files named like this: `control/control/KLF4_MOUSE.PEAKS037825.sissrs.control.mfa`.

These steps listed as separate tasks in Rakefie can be performed using `rake` scripts.
For example to perform `final_collection_summary` task, one should execute `rake final_collection_summary` command.

Some of the tasks work in a different way: they just print other commands that should be run.
One can pipe that commands into bash or (better) to `parallel` (a program from GNU moreutils package).
For example to run 8 threads of logROC AUC calculations in parallel one can invoke the following command:
`rake calculate_occurence_scores_mono | parallel -j 8`

After running all the steps, resulting collection will appear in the folder `final_bundle`.

Supplementary data for figures in paper were prepared using scripts in `data_for_figures.sh`.
