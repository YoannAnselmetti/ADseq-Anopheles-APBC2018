#!/bin/bash
###
###   Goal:
###      Execute the 6 experiments for the validation of ADseq algo with BLASTn alignment
###      	=> 3 species (Aalb, Aara & Adir) x 2 reads sampling (50% & ALL)
###





assemblyTAG="blastn"
stepBegin=1
script="bin/scripts/validation_ADseq/BLASTn/00-run_ALL_script_BLASTn.sh"



#######
### 50pourc
#######
### Anopheles_albimanus
echo ""$script Anopheles_albimanus Aalb 50pourc data/INPUT_DATA/FASTA/SCAFF/Anopheles-albimanus-STECLA_SCAFFOLDS_AalbS1.fa data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_albimanus/50pourc/minia_k75_m3_Aalb_50pourc.contigs.fa $stepBegin AALB $assemblyTAG""
$script Anopheles_albimanus Aalb 50pourc data/INPUT_DATA/FASTA/SCAFF/Anopheles-albimanus-STECLA_SCAFFOLDS_AalbS1.fa data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_albimanus/50pourc/minia_k75_m3_Aalb_50pourc.contigs.fa $stepBegin AALB $assemblyTAG

### Anopheles_arabiensis
echo ""$script Anopheles_arabiensis Aara 50pourc data/INPUT_DATA/FASTA/SCAFF/Anopheles-arabiensis-Dongola_SCAFFOLDS_AaraD1.fa data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_arabiensis/50pourc/minia_k59_m3_Aara_50pourc.contigs.fa $stepBegin AARA $assemblyTAG""
$script Anopheles_arabiensis Aara 50pourc data/INPUT_DATA/FASTA/SCAFF/Anopheles-arabiensis-Dongola_SCAFFOLDS_AaraD1.fa data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_arabiensis/50pourc/minia_k59_m3_Aara_50pourc.contigs.fa $stepBegin AARA $assemblyTAG

### Anopheles_dirus
echo ""$script Anopheles_dirus Adir 50pourc data/INPUT_DATA/FASTA/SCAFF/Anopheles-dirus-WRAIR2_SCAFFOLDS_AdirW1.fa data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_dirus/50pourc/minia_k63_m3_Adir_50pourc.contigs.fa $stepBegin ADIR $assemblyTAG""
$script Anopheles_dirus Adir 50pourc data/INPUT_DATA/FASTA/SCAFF/Anopheles-dirus-WRAIR2_SCAFFOLDS_AdirW1.fa data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_dirus/50pourc/minia_k63_m3_Adir_50pourc.contigs.fa $stepBegin ADIR $assemblyTAG



#######
### ALL
#######
### Anopheles_albimanus
echo ""$script Anopheles_albimanus Aalb ALL data/INPUT_DATA/FASTA/SCAFF/Anopheles-albimanus-STECLA_SCAFFOLDS_AalbS1.fa data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_albimanus/ALL/minia_k83_m3_Aalb_ALL.contigs.fa $stepBegin AALB $assemblyTAG""
$script Anopheles_albimanus Aalb ALL data/INPUT_DATA/FASTA/SCAFF/Anopheles-albimanus-STECLA_SCAFFOLDS_AalbS1.fa data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_albimanus/ALL/minia_k83_m3_Aalb_ALL.contigs.fa $stepBegin AALB $assemblyTAG

### Anopheles_arabiensis
echo ""$script Anopheles_arabiensis Aara ALL data/INPUT_DATA/FASTA/SCAFF/Anopheles-arabiensis-Dongola_SCAFFOLDS_AaraD1.fa data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_arabiensis/ALL/minia_k72_m3_Aara_ALL.contigs.fa $stepBegin AARA $assemblyTAG""
$script Anopheles_arabiensis Aara ALL data/INPUT_DATA/FASTA/SCAFF/Anopheles-arabiensis-Dongola_SCAFFOLDS_AaraD1.fa data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_arabiensis/ALL/minia_k72_m3_Aara_ALL.contigs.fa $stepBegin AARA $assemblyTAG

### Anopheles_dirus
echo ""$script Anopheles_dirus Adir ALL data/INPUT_DATA/FASTA/SCAFF/Anopheles-dirus-WRAIR2_SCAFFOLDS_AdirW1.fa data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_dirus/ALL/minia_k75_m3_Adir_ALL.contigs.fa $stepBegin ADIR $assemblyTAG""
$script Anopheles_dirus Adir ALL data/INPUT_DATA/FASTA/SCAFF/Anopheles-dirus-WRAIR2_SCAFFOLDS_AdirW1.fa data/validation_ADseq/FASTA/SCAFF/minia/Anopheles_dirus/ALL/minia_k75_m3_Adir_ALL.contigs.fa $stepBegin ADIR $assemblyTAG