# Yoann Anselmetti and Cedric Chauve, January 25 2017
# Linearizing the genomes obtained with the validation data sets and new thresholds



DIR=$1
assemblyTAG=$2


ALGORITHM="ADseq ARt-DeCo"
scriptLINEAR="bin/scripts/post_decostar/code/linearize_genomes.py"
for algo in $ALGORITHM; do
	# General commands for the ADseq results
	# In all experiments, the Boltzmann temperature was set at 0.1
	PREF_ADIR_ALL=results/validation_ADseq/$DIR/Adir/ALL/$algo/$assemblyTAG/validation_ADseq_Adir_ALL_Xtopo_Boltz_0.1
	PREF_ADIR_HALF=results/validation_ADseq/$DIR/Adir/50pourc/$algo/$assemblyTAG/validation_ADseq_Adir_50pourc_Xtopo_Boltz_0.1
	PREF_ADIR_ALL1=${PREF_ADIR_ALL}_01
	PREF_ADIR_HALF1=${PREF_ADIR_HALF}_01
	PREF_ADIR_ALL5=${PREF_ADIR_ALL}_05
	PREF_ADIR_HALF5=${PREF_ADIR_HALF}_05
	PREF_ADIR_ALL8=${PREF_ADIR_ALL}_08
	PREF_ADIR_HALF8=${PREF_ADIR_HALF}_08

	echo "linearize_genomes.py ${PREF_ADIR_ALL}.adjacencies.txt 23 0.1 ${PREF_ADIR_ALL1} M1 23 0.01" > ${PREF_ADIR_ALL1}_linearization.log
	python $scriptLINEAR ${PREF_ADIR_ALL}.adjacencies.txt 23 0.1 ${PREF_ADIR_ALL1} M1 23 0.01 >> ${PREF_ADIR_ALL1}_linearization.log 
	echo "linearize_genomes.py ${PREF_ADIR_HALF}.adjacencies.txt 23 0.1 ${PREF_ADIR_HALF1} M1 23 0.01" > ${PREF_ADIR_HALF1}_linearization.log
	python $scriptLINEAR ${PREF_ADIR_HALF}.adjacencies.txt 23 0.1 ${PREF_ADIR_HALF1} M1 23 0.01 >> ${PREF_ADIR_HALF1}_linearization.log

	echo "linearize_genomes.py ${PREF_ADIR_ALL}.adjacencies.txt 23 0.5 ${PREF_ADIR_ALL5} M1 23 0.01" > ${PREF_ADIR_ALL5}_linearization.log
	python $scriptLINEAR ${PREF_ADIR_ALL}.adjacencies.txt 23 0.5 ${PREF_ADIR_ALL5} M1 23 0.01 >> ${PREF_ADIR_ALL5}_linearization.log 
	echo "linearize_genomes.py ${PREF_ADIR_HALF}.adjacencies.txt 23 0.5 ${PREF_ADIR_HALF5} M1 23 0.01" > ${PREF_ADIR_HALF5}_linearization.log
	python $scriptLINEAR ${PREF_ADIR_HALF}.adjacencies.txt 23 0.5 ${PREF_ADIR_HALF5} M1 23 0.01 >> ${PREF_ADIR_HALF5}_linearization.log 

	echo "linearize_genomes.py ${PREF_ADIR_ALL}.adjacencies.txt 23 0.8 ${PREF_ADIR_ALL8} M1 23 0.01" > ${PREF_ADIR_ALL8}_linearization.log
	python $scriptLINEAR ${PREF_ADIR_ALL}.adjacencies.txt 23 0.8 ${PREF_ADIR_ALL8} M1 23 0.01 >> ${PREF_ADIR_ALL8}_linearization.log 
	echo "linearize_genomes.py ${PREF_ADIR_HALF}.adjacencies.txt 23 0.8 ${PREF_ADIR_HALF8} M1 23 0.01" > ${PREF_ADIR_HALF8}_linearization.log
	python $scriptLINEAR ${PREF_ADIR_HALF}.adjacencies.txt 23 0.8 ${PREF_ADIR_HALF8} M1 23 0.01 >> ${PREF_ADIR_HALF8}_linearization.log
done