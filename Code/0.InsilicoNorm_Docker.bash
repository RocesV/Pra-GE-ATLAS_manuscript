#!bin/bash
## trinityrnaseq/trinity docker should be pulled

INPUT=/mnt/f/ATLAS_Pra_2/Pinus_radiata/fastq_rcorrector/input/
OUTPUT=/mnt/f/ATLAS_Pra_2/Pinus_radiata/fastq_rcorrector/output/
DIRS=$(ls ${INPUT})

#  docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /usr/local/bin/util/insilico_read_normalization.pl
# CHECK: careful with the -v volume rewrite all --> first 'pwd' change by /mnt/f/ATLAS_Pra/ and change files to
# relative path of 'pwd'

for j in $DIRS
	do
	echo "Processing  ${j} directory ... | $(date)"
	DIRS2=$(ls ${INPUT}${j})
	for i in $DIRS2
	do
	if [ "${i}" = "PE" ]; then
	echo "Processing paired-end files ..."
	DIRS3=$(ls ${INPUT}${j}/${i})
	for h in $DIRS3
		do
		echo "Processing  ${h} directory ... | $(date)"
		for f in ${INPUT}${j}/${i}/${h}/*_1.cor.fq
			do
			base=${f%_1*}
	  		foward="${base}_1.cor.fq"
	  		reverse="${base}_2.cor.fq"
			if [ "${h}" = "Unstranded" ]; then
			docker run --rm -v /mnt/f/ATLAS_Pra_2/Pinus_radiata/fastq_rcorrector/:/home/ trinityrnaseq/trinityrnaseq /usr/local/bin/util/insilico_read_normalization.pl --seqType fq --JM 60G --max_cov 30 --left /home/input/${j}/${i}/${h}/$(basename ${foward}) --right /home/input/${j}/${i}/${h}/$(basename ${reverse}) --pairs_together --PARALLEL_STATS --CPU 16 --output /home/output/${j}/
			else
			docker run --rm -v /mnt/f/ATLAS_Pra_2/Pinus_radiata/fastq_rcorrector/:/home/ trinityrnaseq/trinityrnaseq /usr/local/bin/util/insilico_read_normalization.pl --seqType fq --JM 60G --max_cov 30 --left /home/input/${j}/${i}/${h}/$(basename ${foward}) --right /home/input/${j}/${i}/${h}/$(basename ${reverse}) --pairs_together --PARALLEL_STATS --CPU 16 --output /home/output/${j}/ --SS_lib_type ${h}
			fi
			done
		done
	else 
	echo "Processing single-end files ..."
	DIRS3=$(ls ${INPUT}${j}/${i})
	for h in $DIRS3
                do
                echo "Processing  ${h} directory ... | $(date)"
		for f in ${INPUT}${j}/${i}/${h}/*
			do
			if [ "${h}" = "Unstranded" ]; then
			docker run --rm -v /mnt/f/ATLAS_Pra/:/home/ trinityrnaseq/trinityrnaseq /usr/local/bin/util/insilico_read_normalization.pl --seqType fq --JM 60G --max_cov 30 --single /home/input/Normalization/${j}/${i}/${h}/$(basename ${f}) --CPU 16 --output /home/output/${j}/
			else
			docker run --rm -v /mnt/f/ATLAS_Pra/:/home/ trinityrnaseq/trinityrnaseq /usr/local/bin/util/insilico_read_normalization.pl --seqType fq --JM 60G --max_cov 30 --single /home/input/Normalization/${j}/${i}/${h}/$(basename ${f}) --CPU 16 --output /home/output/${j}/ --SS_lib_type ${h}
			fi
			done
                done
	fi
	done
	done
