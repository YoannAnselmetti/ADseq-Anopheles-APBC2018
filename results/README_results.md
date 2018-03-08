
Output files format of ADseq algorithm and post-processing steps (adjacencies linearisation, scaffolds assignment and SCJ computations)
===


For all files output by ADseq, the algorithm take as input an argument "output.dir" corresponding to the prefix of output files (cf. an [example](../data/data_DeCoSTAR/decostar/Xtopo_pNJ/DeCoSTAR_Anopheles_Xtopo+scaff.param.txt) config file of DeCoSTAR).
In this example the prefix is: DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1
* Xtopo: topology of the species tree used as input (X topology)
* +scaff: use of sequencing data to compute adjacencies evolutionary histories
* Boltz_0.1: use of Boltzmann distribution with a kT=0.1 that ensures that predictions are highly dominated by optimal and slightly sub-optimal solutions on the parsimonious criterion -> min(#adjacencies breaks and gains). For more details see this [article](https://link.springer.com/chapter/10.1007/978-3-319-12418-6_7).

The ADseq algorithm outputs 5 different files:
*  $(prefix).adjacencies.txt -> each line corresponds to an adjacency (extant or ancestral) outputs by ADseq. Format of the file by  columns:
	1. species ID (correspondence with extant species names is present in files: $(prefix).species.txt and $(prefix).speciesTree.newick)
	2. gene of the adjacency in 5' position -> Extant gene format: $(species_name)$(separator)$(gene_ID), here $(separator) is "@". Ancestral gene format: $(gene_tree_ID)|$(node_ID)
	3. gene of the adjacency in 3' (same formats as above)
	4. strand of the 5' gene (+: forward and -: reverse)
	5. strand of the 3' gene
	6. a priori support. Here, if support is equal to 1 then adjacency was present in initial assembly, if >0 and <=1 then adjacency is supported by paired sequencing data, and if equal to 0 the adjacency was not present and not supported by sequencing data)
	7. a posteriori support corresponds to the ratio of sampled evolutionary scenarios where the adjacency is present on the total number of sampled evolutionary scenarios.
*  $(prefix).genes.txt -> each line corresponds to a gene corresponding to a speciation event or a leaf (not duplicated and loss genes) considered by ADseq (i.e. present in reconciled gene trees). Format of the file by columns:
	1. species ID
	2. gene ID -> Extant gene format: $(species_name)$(separator)$(gene_ID), here $(separator) is "@". Ancestral gene format: $(gene_tree_ID)|$(node_ID)
	3. columns 3 to the last ones corresponds to child of the gene corresponding to speciation or leaf node. If one of the child is a duplicated then take the children of this duplication.
*  $(prefix).reconciliations.newick: reconciled gene trees considered by ADseq to computed adjacencies evolutionary histories. Node format: $(gene_ID)|$(evolutionary_event)|$(species_ID)|$(time_slice). Here, $(time_slice) is not valid cause we don't take into account horizontal gene transfer.
*  $(prefix).species.txt: association between species ID and species name
*  $(prefix).speciesTree.newick: species tree with species ID (and species name for extant ones)


After DeCoSTAR, adjacencies have to be linearised due to independent computation of each adjacency history resulting in syntenic conflicts between inferred adjacencies. These steps results in 11 output files with the same $(prefix) used for ADseq output files, with linearisation with M1 algorithm (see this [script](../bin/scripts/post_decostar/code/linearize_genomes.py) for more details) with a adjacencies with a posteriori support superior or equal to 0.1:  
* $(prefix)\_0.1_M1_disc: list of adjacencies that have been discarded during linearisation step. Last column "filtering_stage" has two possible values: "1" corresponds to the adjacencies that have been discarded due to a posterior support under fixed threshold (here, 0.1). "MWM" corresponds to the adjacencies filtered during the Maximum-Weight-Matching process.
* $(prefix)\_0.1_M1_kept: list of adjacencies that have been kept after linearisation process. Same format as previous file except that last column corresponds to the scaffold ID inferred by ADseq.
* $(prefix)\_0.1_M1_new_extant_adjacencies: list of extant adjacencies that have been inferred by ADseq and passed linearisation process
* $(prefix)\_0.1_M1_obs_scaffolds: list of genes with strand orientation, species and scaffolds ID they belong in initial genome assembly
* $(prefix)\_0.1_M1_scaffolds: list of genes with strand orientation, species and scaffolds ID they belong after ADseq scaffolding and linearisation process
* $(prefix)\_0.1_M1_scaffolds_assignment: list of scaffolds/gene assignment to the Anopheles chromosomes (X, 2L, 2R, 3L, 3R)
* $(prefix)\_0.1_M1_scaffolds_log: log file of adjacencies linearisation
* $(prefix)\_0.1_M1_scj: list of SCJ for each internal branches of the species tree. 2 SCJ 
* $(prefix)\_0.1_M1_scj_log: log file of SCJ computations
* $(prefix)\_0.1_M1_scj_stats: file summarizing statistics on SCJ computations
* $(prefix)\_0.1_orthogroups -> 1 line per pair of ancestor-descendant gene, branches are only the pre-extant or pre-speciation branches. File format by columns:
	1. ancestor species ID
	2. descendant species ID
	3. ancestor gene ID
	4. descendant gene ID
	5. gene tree ID
