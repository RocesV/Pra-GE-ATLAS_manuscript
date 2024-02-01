#!bin/bash
TOOLS=/mnt/g/ATLAS_Pra/tools/Trimmomatic-0.39/
INPUT=/mnt/g/ATLAS_Pra/ATLAS/fastq_sortmerna/
OUTPUT=/mnt/g/ATLAS_Pra/ATLAS/fastq_trimmomatic/
SUMMARY=/mnt/g/ATLAS_Pra/ATLAS/summary_trimmomatic/
DIRS=$(ls ${INPUT})


for j in $DIRS
	do
	echo "Processing ${j} directory ... | $(date)"
	if [ "${j}" = "PE" ]; then
	echo "Processing paired-end files ..."
	 for f in ${INPUT}${j}/*_1.fastq
	  do
	  base=${f%_1*}
	  foward="${base}_1.fastq"
	  reverse="${base}_2.fastq"
      echo "Trimming $base files ... | $(date)"
	  echo "Analysis phase"
	  java -jar ${TOOLS}trimmomatic-0.39.jar PE -threads 16 -phred33 -summary ${SUMMARY}$(basename ${base})_analysis_summary.txt ${foward} ${reverse} ${OUTPUT}PE/analysis/Paired_$(basename ${foward}) ${OUTPUT}PE/analysis/Unpaired_$(basename ${foward}) ${OUTPUT}PE/analysis/Paired_$(basename ${reverse}) ${OUTPUT}PE/analysis/Unpaired_$(basename ${reverse}) ILLUMINACLIP:${TOOLS}adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:48 
	  echo "Assembly phase"
	  java -jar ${TOOLS}trimmomatic-0.39.jar PE -threads 16 -phred33 -summary ${SUMMARY}$(basename ${base})_assembly_summary.txt ${foward} ${reverse} ${OUTPUT}PE/assembly/Paired_$(basename ${foward}) ${OUTPUT}PE/assembly/Unpaired_$(basename ${foward}) ${OUTPUT}PE/assembly/Paired_$(basename ${reverse}) ${OUTPUT}PE/assembly/Unpaired_$(basename ${reverse}) ILLUMINACLIP:${TOOLS}adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:5:5 MINLEN:48 	
	  done
	else
	echo "Processing single-end files ..."
	 for f in ${INPUT}${j}/*
	  do
	  echo "Trimming $f files ... | $(date)"
	  echo "Analysis phase"
	  java -jar ${TOOLS}trimmomatic-0.39.jar SE -threads 16 -phred33 -summary ${SUMMARY}$(basename ${f})_analysis_summary.txt ${f} ${OUTPUT}SE/analysis/$(basename ${f}) ILLUMINACLIP:${TOOLS}adapters/TruSeq2-SE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:48
	  echo "Assembly phase"
	  java -jar ${TOOLS}trimmomatic-0.39.jar SE -threads 16 -phred33 -summary ${SUMMARY}$(basename ${f})_assembly_summary.txt ${f} ${OUTPUT}SE/assembly/$(basename ${f}) ILLUMINACLIP:${TOOLS}adapters/TruSeq2-SE.fa:2:30:10 SLIDINGWINDOW:5:5 MINLEN:48
	  done
	fi 
	done 
