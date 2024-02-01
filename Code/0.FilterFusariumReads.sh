#!bin/bash

INPUT=/mnt/e/ATLAS_Pra/ATLAS/Pinus_radiata/z.fastq_fusariumgenome_filtered/rcorrector/
OUTPUT=/mnt/e/ATLAS_Pra/ATLAS/Pinus_radiata/z.fastq_fusariumgenome_filtered/rcorrector_filtered/
INDEX=/mnt/e/ATLAS_Pra/ATLAS/Pinus_radiata/z.starfiles/index/FSP34
DIRS=$(ls ${INPUT})

for j in $DIRS
	do
	echo "Processing ${j} directory ... | $(date)"
	if [ "${j}" = "PE" ]; then
	echo "Processing paired-end files ..."
	 for f in ${INPUT}${j}/*_1.cor.fq
	  do
	  base=${f%_1*}
	  foward="${base}_1.cor.fq"
	  reverse="${base}_2.cor.fq"
          echo "Fusarium genome-based filtering $base files ... | $(date)"
          bowtie2 -q --end-to-end --un-conc ${OUTPUT}/$(basename ${base}) -p 16 -x ${INDEX} -1 ${foward} -2 ${reverse}
	  done
	else
	echo "Processing single-end files ..."
	 for f in ${INPUT}${j}/*
	  do
	  echo "Fusarium genome-based filtering $f files ... | $(date)"
          bowtie2 -q --end-to-end --un ${OUTPUT}/$(basename ${f}) -p 16 --quiet -x ${INDEX} -U ${f}
	  done
	fi
	done
