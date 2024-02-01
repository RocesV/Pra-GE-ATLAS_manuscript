#!bin/bash
TOOLS=/mnt/e/ATLAS_Pra/tools/rcorrector/
INPUT=/mnt/e/ATLAS_Pra/ATLAS/fastq_trimmomatic/Assembly/
OUTPUT=/mnt/d/0.ATLAS_Pra/fastq_rcorrector/
DIRS=$(ls ${INPUT})

for j in $DIRS
	do
	echo "Processing ${j} directory ... | $(date)"
	if [ "${j}" = "Paired" ]; then
	echo "Processing paired-end files ..."
	 for f in ${INPUT}${j}/*_1.fastq
	  do
	  base=${f%_1*}
	  foward="${base}_1.fastq"
	  reverse="${base}_2.fastq"
      echo "RCorrecting $base files ... | $(date)"
	  echo "Lets go young boy"
	  perl ${TOOLS}run_rcorrector.pl -1 ${foward} -2 ${reverse} -t 18 -od ${OUTPUT}Paired/
	  done
	else
	echo "Processing single-end files ..."
	 for f in ${INPUT}${j}/*
	  do
	  echo "RCorrecting $f files ... | $(date)"
	  echo "Lets go young boy"
	  perl ${TOOLS}run_rcorrector.pl -s ${f} -t 18 -od ${OUTPUT}Single/
	  done
	fi 
	done 
