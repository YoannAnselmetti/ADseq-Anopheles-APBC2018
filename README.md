Phylogenetic signal from rearrangements in 18 Anopheles species by joint scaffolding extant and ancestral genomes
=====

# Introduction

This repository contains all data and results produced for the paper "Phylogenetic signal from rearrangements in 18 Anopheles species by joint scaffolding extant and ancestral genomes" (accepted to [APBC2018](http://apbc2018.bio.keio.ac.jp/) and to appear in BMC Genomics). This paper presents a new method, called ADseq, that allows to improve scaffolding of extant genomes and jointly reconstructs gene order in ancestral genomes. It is mainly focused on the improvement of scaffolding of 18 Anopheles extant genomes.

First download or clone the GitHub repository with the green button "Clone or download" or with the command line:

```
git clone https://github.com/YoannAnselmetti/ADseq-Anopheles-APBC2018.git
```


# Repository structure and content

The scheme below represents the global structure of this repository: 
.
├── bin/
├── data/
├── doc/18Anopheles_sequencing_data.ods
├── figures/
├── README.md
└── results/

In the following, we describe the different directories present in this GitHub repository. 


## "bin/" directory
bin
├── scripts
│   ├── before_decostar
│   │   ├── code/
│   │   ├── run_all_distrib_gene_trees.sh
│   │   └── run_all_plot_hist_scaff_BESST.sh
│   ├── clean_RAW_DATA
│   │   ├── clean_INPUT_DATA.sh
│   │   └── code/
│   ├── post_decostar
│   │   ├── code/
│   │   ├── compute_stats_graph/
│   │   │   ├── code/
│   │   │   ├── run_all_ETE.sh
│   │   │   ├── run_all_newADJfiles.sh
│   │   │   ├── run_all_R_plots.sh
│   │   │   └── run_all_scatterplots.sh
│   │   ├── linearize_generate_stats_decostar.sh
│   │   └── others.zip
│   ├── run_decostar.sh
│   └── validation_ADseq/
│       ├── 01-run_ALL_exp_minia.sh
│       ├── 02-run_ALL_exp_BLASTn.sh
│       ├── 03-run_decostar_validation_and_linearization.sh
│       ├── 04-stats_graphics_validation.sh
│       ├── BLASTn/
│       ├── decostar/
│       ├── minia/
│       └── stats/
└── software_libraries/

The directory **bin/** contains the softwares and scripts that have been used and developed to process raw input data of the 18 Anopheles dataset (mainly produced by [Neafsey et al., 2015](http://science.sciencemag.org/content/347/6217/1258522.long), located in directory **data/INPUT_DATA** and described in section **data/** directory to produce results present in the paper.

* The directory **bin/software_libraries** contains the libraries and the executive files of softwares (such as DeCoSTAR, BESST, minia, SAMtools, ...) that have been used in the experiments of our paper. To reproduce the results present in this repository it is important to uncompressed zipped directories present in directory "bin/software_libraries".
* The directory **bin/scripts** contains all scripts use to produce results of  our paper. The different scripts have to be executed in a certain order to reproduce the results form input raw data files present in directory **data/INPUT_DATA** (cf section **"data/ directory**):
    1. Execute the script **bin/scripts/clean_RAW_DATA/clean_INPUT_DATA.sh** that uses two scripts present in directory **bin/scripts/clean_RAW_DATA/code** to clean raw input data and format them in standard format use by the pipeline that will produce input data for the DeCoSTAR software.
    
    2. Use the pipeline to produce input data for DeCoSTAR (available in the GitHub repository [DeCoSTAR_pipeline](https://github.com/YoannAnselmetti/DeCoSTAR_pipeline)) that will generate, from the files present in the **data/INPUT_DATA/** directory, the input data for the DeCoSTAR software that will be stored in the directory **data/**. This pipeline contains all the steps illustarte in Figure 6 of our paper. However some data (like data sequencing files (FASTQ format)) are not available due to excessive size of the files. For now the pipeline doesn't allow to execute the pipeline to infer gene trees with RAxML and profileNJ and the pipeline to produce scaffolding adjacencies with BESST (They have to be implemented in snakemake to be executable on a cluster with SGE architecture). Shortly, the **DeCoSTAR_pipeline** GitHub repository will get further developments to allow a full run of DeCoSTAR (input data production for DeCoSTAR, DeCoSTAR execution, linearization of adjacencies predictions and computation of graph/figures statistics). 
    
    3. The directory **bin/scripts/validation_ADseq/** contains scripts developed to validate scaffolding ability of the algorithm ADseq (cf subsection "Validation of the ADseq algorithm for extant scaffolding" of the paper). The main idea is to fragment the genome of 3 species located with different characteristics in the species tree:
    * A. albimanus: outgroup position
    * A. arabiensis: deep in the species tree with several close related species
    * A. dirus: deep in the species tree with few close related species
    
        by producing a more conservative genome assembly with the assembly tool minia with 2 different reads sampling (50% and 100%).
This validation process is divided in 4 parts:
        1. produce new genome assembly with minia for the 3 species and the 2 reads sampling and scaffold minia  adjacencies with BESST (with script **bin/scripts/validation_ADseq/01-run_ALL_exp_minia.sh** using scripts present in directory **bin/scripts/validation_ADseq/minia**). This script corresponds to the steps 1/, 2/ and 3/ of Figure 12 of our paper. These steps are not reproducible from this repository due to missing data, corresponding to the sequencing data in the FASTQ format that are too heavy to be stored on a GitHub repository. 
        2. map minia contigs on reference genome with BLASTn for the 6 experiments to transfer gene annotation from the reference genome to the minia contigs (with script **bin/scripts/validation_ADseq/02-run_ALL_exp_BLASTn.sh** using scripts present in directory **bin/scripts/validation_ADseq/BLASTn**). This script corresponds to the steps 4/ and 5/ f Figure 12. However, the script doesn't contain a step to execute BESST on BLASTn scaffolds to compute scaffolding adjacencies that has to be done after the script **bin/scripts/validation_ADseq/BLASTn/02-filter_BLASTn_results** and before the script **bin/scripts/validation_ADseq/BLASTn/06a-create_scaff_adj_prefile_BLASTn.py**
        3. execute DeCoSTAR and linearize adjacencies predictions (with script **bin/scripts/validation_ADseq/03-run_decostar_validation_and_linearization.sh** using scripts present in directory **bin/scripts/validation_ADseq/decostar**)
        4. compare adjacencies predictions of ARt-DeCo, ADseq and BESST to the reference genome (precision and recall statistics) and compare results between the 3 methods: Venn diagrams (with script **bin/scripts/validation_ADseq/04-stats_graphics_validation.sh** using scripts present in directory **bin/scripts/validation_ADseq/stats**) 

    4. The directory **bin/scripts/post_decostar** contains the script **linearize_generate_stats_decostar.sh** that linearizes adjacencies predicted by DeCoSTAR and compute statistics on genome rearrangements and scaffolding inferred and produced by DeCoSTAR with scripts present in directory **bin/scripts/post_decostar/code**. It contains also the directory **compute_stats_graph/**  that allows to compute all graphs and figures on the statistics of the results of DeCoSTAR present in the paper. 



## "data/" directory
data
├── data_DeCoSTAR
├── DATA_SEQ
├── FASTA
├── GENE_TREES
├── GFF_to_GENE_files
├── INPUT_DATA
└── validation_ADseq

The **data/ ** directory contains all files used to produce input data to apply the DeCoSTAR software on the 18 Anopheles dataset.


### "data/INPUT_DATA/" directory
data/INPUT_DATA/
├── 18Anopheles_species
├── Anopheles_species_tree_WG_topology.nwk
├── Anopheles_species_tree_X_topology.nwk
├── FASTA
├── GFF
├── name_geneID_18Anopheles
├── OG_CDS_newtrees.tar.gz
├── ORIGINAL_GFF
└── unrooted_raw_trees.nwk

The **data/INPUT_DATA** directory contains all input data available for the 18 Anopheles dataset mainly produced by [Neafsey et al., 2015](http://science.sciencemag.org/content/347/6217/1258522.long). Among the 18 Anopheles genomes, 4 have no paired sequencing data available (annotate with "X") to produce scaffolding adjacencies for ADseq algorithm. For more details on sequencing data available see table file **doc/18Anopheles_sequencing_data.ods**:

* Anopheles_albimanus
* Anopheles_arabiensis
* Anopheles_atroparvus
* Anopheles_christyi
* Anopheles_culicifacies
* Anopheles_darlingi        X
* Anopheles_dirus
* Anopheles_epiroticus
* Anopheles_farauti
* Anopheles_funestus
* Anopheles_gambiae     X
* Anopheles_maculatus
* Anopheles_melas
* Anopheles_merus
* Anopheles_minimus
* Anopheles_quadriannulatus X
* Anopheles_sinensis
* Anopheles_stephensi       X

The gene annotation files in GFF format (directory **data/INPUT_DATA/ORIGINAL_GFF**), the unrooted gene trees in newick format (directory **data/INPUT_DATA/OG_CDS_newtrees**) and the CDS (directory **data/INPUT_DATA/FASTA/CDS**) for the 18 Anopheles species have been obtained from Robert M. Waterhouse. The 18 Anopheles reference genome assemblies have been obtained from [VectorBase](https://www.vectorbase.org/downloads) (directory **data/INPUT_DATA/FASTA/SCAFF**). Original GFF files have been processed with script **bin/scripts/clean_RAW_DATA/clean_INPUT_DATA.sh** to format them to the standard GFF format file and have been stored in the directory **data/INPUT_DATA/GFF**.  
There are also 4 files in **data/INPUT_DATA**:

* **18anopheles_species**: handmade file associating species name with number of chromosome expected
* **Anopheles_species_tree_X_topology.nwk**: 18 Anopheles species tree with the X chromosome topology (X)
* **Anopheles_species_tree_WG_topology.nwk**: 18 Anopheles species tree with the Whole Genome topology (WG)
* **name_geneID_18Anopheles**: handmade file associating species name with species ID
* **unrooted_raw_trees.nwk**; contains the gene trees after cleaning of gene trees available in the directory **data/INPUT_DATA/OG_CDS_newtrees** with the script **bin/scripts/clean_RAW_DATA/clean_INPUT_DATA.sh**



### "data/GFF_to_GENE_files/" directory
data/GFF_to_GENE_files
├──filtered_GENE
├── GENE
├── GRAPH_GFF
├── sorted_GENE
├── sorted_GFF
└── with_filter


The **data/GFF_to_GENE_files/** directory contains intermediate file to transform GFF files in GENE file.
Initial GFF files located in directory **data/INPUT_DATA/GFF** are sorted by gene and exon positions in directory **data/GFF_to_GENE_files/sorted_GFF**. Then, sorted GFF files are transformed in GENE files in directory **data/GFF_to_GENE_files/GENE_file**. They are then sorted in directory **data/GFF_to_GENE_files/sorted_GENE**, sorted GENE files are restricted to genes present in gene trees considered and stored in directory **data/GFF_to_GENE_files/filtered_GENE**. Genes positions are analyzed to produce 3 files are stored in the directory **data/GFF_to_GENE_files/with_filter**:

* **ALL_species_Inclusion_file**: file containing included genes
*  **ALL_species_Overlap_file**: file containing overlapping genes
*   **ALL_species_GENE_file**: file containing non included genes of the 18 Anopheles species

Finally, a last file is produce to add gene family ID to GENE file (**data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file_with_GF**)



### "data/FASTA" directory
data/FASTA
├── MSA/CDS
│    ├── Gblocks
│    └── MUSCLE
└──  GF_FASTA


The **FASTA** directory contains data to infer gene trees with the gene tree pipeline inference described in our paper with the refinement gene tree tool profileNJ.

The directory **data/FASTA/GT_FASTA** contains 1 FASTA files / gene family/tree containing the CDS of the genes belonging to the gene family/tree.
The **data/FASTA/MSA/CDS** contains two directories:
    * **MUSCLE** which contains Multiple Sequence Alignment (MSA) files for each gene family produced with the MSA tool MUSCLE. 
    * **Gblocks** which contains selected blocks of multiple sequence alignment produced with the tool Gblocks.



### "data/GENE_TREES" directory
data/GENE_TREES
├── CDS
│   └── bootstrap_support
│       ├── profileNJ
│       │   ├── UNROOTED_GENE_TREES.tar.gz
│       │   ├── WG_topo.tar.gz
│       │   └── X_topo.tar.gz
│       └── RAxML
├── trees_DeCoSTAR_WGtopo.nwk
├── trees_DeCoSTAR_Xtopo.nwk
└── unrooted_trees_filtered.nwk


The **data/GENE_TREES** directory contains the gene trees files produced for the 18 Anopheles dataset.
The directory contains 6 files:
    * **trees_DeCoSTAR_WGtopo.nwk**: gene trees used as input of DeCoSTAR with the species tree topology WG
    * **trees_DeCoSTAR_Xtopo.nwk**: gene trees used as input of DeCoSTAR with the species tree topology 
It contains also 2 directories:
    * **CDS/bootstrap_support/RAxML** contains gene trees that have been produced with the maximum likelihood inference gene tree tool RAxML.
    * **CDS/bootstrap_support/profileNJ** contains gene trees that have been produced with the refinement gene tree tool profileNJ.


### "data/DATA_SEQ" directory
data/DATA_SEQ
├── orientation_libraries
└── SCAFFOLDING/BESST-2.2.6
    ├── Bowtie2_ALL/TRIMMOMATIC3/ALL/
    └── Bowtie2_k50/TRIMMOMATIC3/blastn
        ├── 50pourc
        └── ALL

The **data/DATA_SEQ** directory contains scaffolding results obtained with the scaffolding tool BESST. The directory contain one file **orientation_libraries** that contains the information on the orientation of paired reads for the different SRX (cf file **doc/18Anopheles_sequencing_data.ods**) that are necessary to execute BESST. The FASTQ and BAM files are not present in this repository cause they take too much place (several To of space memory).
The directory **data/DATA_SEQ/SCAFFOLDING/BESST-2.2.6** contains two directories.
    * **Bowtie2_ALL/TRIMMOMATIC3/ALL/** contains scaffolding results of BESST on the 14 Anopheles reference genomes for which sequencing data are available.
    * **Bowtie_k50/TRIMMOMATIC3/blastn** contains scaffolding results of BESST on the de novo genome assemblies produced with minia to validate the ability of ADseq to scaffold genomes with 50% of the reads (directory **50pourc**) and with all reads (directory **ALL**).


### "data/data_DeCoSTAR" directory
data/data_DeCoSTAR/
├── CTG_file
├── decostar
│   ├── adjacencies.txt
│   ├── WGtopo+scaff
│   ├── Xtopo_pNJ
│   └── Xtopo_RAW
├── GENE_file
├── scaff_BESST_ALL_3_TRIMMOMATIC3
└── scaff_BESST_DeCoSTAR

The **data/data_DeCoSTAR/** directory contains input data files to execute DeCoSTAR on the 18 Anopheles dataset:
* The **CTG_file** contains informations on all contigs/scaffolds considered as input of DeCoSTAR for the 18 Anopheles species.
* The **GENE_file** contains informations on all genes considered as input of DeCoSTAR for the 18 Anopheles species.
* The **scaff_BESST_ALL_3_TRIMMOMATIC3** contains the scaffolding adjacencies between contigs/scaffolds of reference genome assemblies computed by BESST with 2 link scores for the 18 Anopheles species.
* The **scaff_BESST_ALL_DeCoSTAR** contains the scaffolding adjacencies between contigs/scaffolds of computed by BESST with 2 link scores between contigs/scaffolds considered as input of DeCoSTAR.
* The **data/data_DeCoSTAR/decostar/** directory contains files use in the parameter files to execute DeCoSTAR: The **adjacencies.txt** file contains all adjacencies considered by DeCoSTAR obtained from the files **data/data_DeCoSTAR/GENE_file** and **data/data_DeCoSTAR/scaff_BESST_DeCoSTAR**. The 3 directories: **WGtopo+scaff**, **Xtopo_pNJ** and **Xtopo_RAW** contains parameters files and gene trees to apply DeCoSTAR for the 4 different experiments described in our paper.



### "data/validation_ADseq" directory
data/validation_ADseq/
├── BLASTn
│   ├── Anopheles_albimanus/
│   ├── Anopheles_arabiensis/
│   └── Anopheles_dirus/
├── DeCoSTAR/BLASTn/
│   ├── Anopheles_albimanus/
│   ├── Anopheles_arabiensis/
│   └── Anopheles_dirus/
├── FASTA/SCAFF
│   ├── blastn/
│   └── minia/
└── kmergenie
    ├── Anopheles_albimanus/
    ├── Anopheles_arabiensis/
    └── Anopheles_dirus/
    
The **data/validation_ADseq/** directory contains input data files to execute DeCoSTAR with the ADseq and the ARt-DeCo algorithms to validate the ability of the ADseq algorithm to scaffold genomes. The directory is composed of 4 directories:
* The **kmergenie/** directory contains results of the tool Kmergenie to define the best kmer size to assemble reads sampling with minia genome assembly tool.
* The **FASTA/SCAFF** directory contains genome assemblies in the FASTA file format for the 3 species selected for the validation steps obtain with minia (**minia/** directory) and after BLASTn gene annotation (**blastn/** directory).
* The **BLASTn/** directory contains the alignments of contigs produced with minia and maaped on  reference genome assemblies with BLASTn to map gene of reference genome on minia contigs for the 6 experiments (3 species and 2 reads sampling)  
*  The **DeCoSTAR/BLASTn/** directory contains the data to apply DeCoSTAR with the ADseq and ARt-DeCo algorithms on the 6 genome fragmentation experiments.  



## "results" directory
results
├── decostar
│   ├── WGtopo+scaff
│   ├── Xtopo_RAW
│   ├── Xtopo-scaff
│   └── Xtopo+scaff
└── validation_ADseq
    └── spi_20
        ├── Aalb
        ├── Aara
        ├── Adir
        └── stats/blastn

The **results/** directory contains results produced after DeCoSTAR execution.

The **results/decostar** directory contains the results of DeCoSTAR with the algorithm ADseq on the 18 Anopheles dataset for the 4 conditions described in the paper:
* **Xtopo+scaff**: X species tree topology with scaffolding adjacencies and profileNJ trees
* **WGtopo+scaff**: WG species tree topology with scaffolding adjacencies and profileNJ trees
* **Xtopo-scaff**: X species tree topology without scaffolding adjacencies and profileNJ trees
* **Xtopo_RAW**: X species tree topology with scaffolding adjacencies and original gene trees

The **results/validation_ADseq** directory contains the results of DeCoSTAR with the algorithms ADseq and ARt-DeCo for the 6 experiments with a fragmented genome for the validation of the scaffolding ability of ADseq algorithm.
The **spi_20/stats/blastn/** directory contains files on precision and recall statistics for the 6 fragmented genomes experiments for the 3 methods coppared (ADseq, ARt-DeCo and BESST)    and files to produce Venn diagrams to compare scaffolding results between the 3 methods. 



## "figures" directory

figures
├── besst_score
├── distrib_gene_trees
├── ETE_species_trees
├── precision_recall
├── R_plots
├── scatterplot
└── Venn_diagram

The directory **figures/** contains statistics graphs and figures present in the paper and produce from files of **data/** directory and **results/** directory.





# Software used in this study
* [DeCoSTAR](http://pbil.univ-lyon1.fr/software/DeCoSTAR/) - DeCoSTAR software (containing ARt-DeCo and ADseq algorithms). [GitHub repository](https://github.com/WandrilleD/DeCoSTAR)
* [BESST](https://github.com/ksahlin/BESST) - BESST scaffolding tool
* [profileNJ](https://github.com/maclandrol/profileNJ) - profileNJ refinement gene tree tool
* [Samtools](http://samtools.sourceforge.net/)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [minia](http://minia.genouest.org/)
* [BLASTn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Nucleotides&PROGRAM=blastn&BLAST_PROGRAMS=blastn&PAGETYPE=BlastSearch&DATABASE=refseq_rna&DESCRIPTIONS=100&EQ_TEXT=arabidopsis[orgn]&QUERY=8033)
* [GBlocks](http://molevol.cmima.csic.es/castresana/Gblocks.html)
* [Kmergenie](http://kmergenie.bx.psu.edu/)
* [Muscle](https://www.drive5.com/muscle/)
* [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html)
* [SRAtoolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)



# Authors
* **Yoann ANSELMETTI**
* **Wandrille DUCHEMIN**
* **Éric TANNIER**
* **Cedric CHAUVE**
* **Sèverine BÉRARD**
