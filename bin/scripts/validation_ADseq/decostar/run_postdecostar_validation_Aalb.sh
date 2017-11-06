# Yoann Anselmetti and Cedric Chauve, January 25 2017
# Linearizing the genomes obtained with the validation data sets and new thresholds



DIR=$1
assemblyTAG=$2


ALGORITHM="ADseq ARt-DeCo"
scriptLINEAR="bin/scripts/post_decostar/code/linearize_genomes.py"
for algo in $ALGORITHM; do
	# General commands for the ADseq results
	# In all experiments, the Boltzmann temperature was set at 0.1
	PREF_AALB_ALL=results/validation_ADseq/$DIR/Aalb/ALL/$algo/$assemblyTAG/validation_ADseq_Aalb_ALL_Xtopo_Boltz_0.1
	PREF_AALB_HALF=results/validation_ADseq/$DIR/Aalb/50pourc/$algo/$assemblyTAG/validation_ADseq_Aalb_50pourc_Xtopo_Boltz_0.1
	PREF_AALB_ALL1=${PREF_AALB_ALL}_01
	PREF_AALB_HALF1=${PREF_AALB_HALF}_01
	PREF_AALB_ALL5=${PREF_AALB_ALL}_05
	PREF_AALB_HALF5=${PREF_AALB_HALF}_05
	PREF_AALB_ALL8=${PREF_AALB_ALL}_08
	PREF_AALB_HALF8=${PREF_AALB_HALF}_08

	echo "linearize_genomes.py ${PREF_AALB_ALL}.adjacencies.txt 31 0.1 ${PREF_AALB_ALL1} M1 31 0.01" > ${PREF_AALB_ALL1}_linearization.log
	python $scriptLINEAR ${PREF_AALB_ALL}.adjacencies.txt 31 0.1 ${PREF_AALB_ALL1} M1 31 0.01 >> ${PREF_AALB_ALL1}_linearization.log 
	echo "linearize_genomes.py ${PREF_AALB_HALF}.adjacencies.txt 31 0.1 ${PREF_AALB_HALF1} M1 31 0.01" > ${PREF_AALB_HALF1}_linearization.log
	python $scriptLINEAR ${PREF_AALB_HALF}.adjacencies.txt 31 0.1 ${PREF_AALB_HALF1} M1 31 0.01 >> ${PREF_AALB_HALF1}_linearization.log

	echo "linearize_genomes.py ${PREF_AALB_ALL}.adjacencies.txt 31 0.5 ${PREF_AALB_ALL5} M1 31 0.01" > ${PREF_AALB_ALL5}_linearization.log
	python $scriptLINEAR ${PREF_AALB_ALL}.adjacencies.txt 31 0.5 ${PREF_AALB_ALL5} M1 31 0.01 >> ${PREF_AALB_ALL5}_linearization.log 
	echo "linearize_genomes.py ${PREF_AALB_HALF}.adjacencies.txt 31 0.5 ${PREF_AALB_HALF5} M1 31 0.01" > ${PREF_AALB_HALF5}_linearization.log
	python $scriptLINEAR ${PREF_AALB_HALF}.adjacencies.txt 31 0.5 ${PREF_AALB_HALF5} M1 31 0.01 >> ${PREF_AALB_HALF5}_linearization.log 

	echo "linearize_genomes.py ${PREF_AALB_ALL}.adjacencies.txt 31 0.8 ${PREF_AALB_ALL8} M1 31 0.01" > ${PREF_AALB_ALL8}_linearization.log
	python $scriptLINEAR ${PREF_AALB_ALL}.adjacencies.txt 31 0.8 ${PREF_AALB_ALL8} M1 31 0.01 >> ${PREF_AALB_ALL8}_linearization.log 
	echo "linearize_genomes.py ${PREF_AALB_HALF}.adjacencies.txt 31 0.8 ${PREF_AALB_HALF8} M1 31 0.01" > ${PREF_AALB_HALF8}_linearization.log
	python $scriptLINEAR ${PREF_AALB_HALF}.adjacencies.txt 31 0.8 ${PREF_AALB_HALF8} M1 31 0.01 >> ${PREF_AALB_HALF8}_linearization.log
done