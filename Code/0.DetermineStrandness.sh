#!bin/bash
## From conda env RepsMask HomePCWork

BED=/mnt/d/Splicing_codan/Pra_res/annotation.bed
INDEX=/mnt/d/1.strandness_wd/index/EPIOMIDES
OUTPUT=/mnt/d/1.strandness_wd/output/
INPUT=/mnt/d/0.ATLAS_Pra/Pinus_radiata/0.Whole_Assembly/1.fastq_assembly/strandness_todo/
GROUP=$(ls ${INPUT})

for g in $GROUP
	do
	echo "Processing $g files ... | $(date)"
	LIB=$(ls ${INPUT}${g})
	for l in $LIB
		do
		echo "		Processing $g/$l files ... | $(date)"
		if [ "${l}" == "Paired" ]; then
		for f in ${INPUT}/${g}/${l}/*_1.cor.fq
			do
			base=${f%_1*}
			foward="${base}_1.cor.fq"
	  		reverse="${base}_2.cor.fq"
			echo " 			Processing $base files ... | $(date)"
			# script for paired
			bowtie2 -q --end-to-end -p 16 -x ${INDEX} -S ${OUTPUT}$(basename ${base}).sam -1 ${foward} -2 ${reverse}
			samtools view -bS ${OUTPUT}$(basename ${base}).sam > ${OUTPUT}$(basename ${base}).bam
			samtools sort ${OUTPUT}$(basename ${base}).bam -o ${OUTPUT}$(basename ${base}).bam
			infer_experiment.py -r ${BED} -i ${OUTPUT}$(basename ${base}).bam > ${OUTPUT}$(basename ${base}).txt
			rm ${OUTPUT}*.bam
			rm ${OUTPUT}*.sam
			done
		elif [ "${l}" == "Single" ]; then
		for f in ${INPUT}/${g}/${l}/*
			do
			echo " 			Processing  $f files ... | $(date)"
			# script for single
			bowtie2 -q --end-to-end -p 16 -x ${INDEX} -S ${OUTPUT}$(basename ${f}).sam -U ${f}
			samtools view -bS ${OUTPUT}$(basename ${f}).sam > ${OUTPUT}$(basename ${f}).bam
			samtools sort ${OUTPUT}$(basename ${f}).bam -o ${OUTPUT}$(basename ${f}).bam
			infer_experiment.py -r ${BED} -i ${OUTPUT}$(basename ${f}).bam > ${OUTPUT}$(basename ${f}).txt 
			rm ${OUTPUT}*.bam
                        rm ${OUTPUT}*.sam
			done
		else
			# no script for unpaired
			echo "		CONTINUE"
			continue
		fi
		done
	done
