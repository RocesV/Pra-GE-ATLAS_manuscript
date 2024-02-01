#!bin/bash

INDEX=/mnt/g/Pra-GE-ATLAS/ATLAS/Pinus_radiata/1.Assemblies/3.Corset/Corset_okayaltokayReduced/0.mapping_bams/salmon/index/Pra-GE-ATLAS_Pra_alloksReduced
OUTPUT=/mnt/g/Pra-GE-ATLAS/ATLAS/Pinus_radiata/1.Assemblies/3.Corset/Corset_okayaltokayReduced/0.mapping_bams/salmon/mapping/
INPUT=/mnt/g/Pra-GE-ATLAS/ATLAS/Pinus_radiata/0.Whole_Analysis/salmon_todo/RF/
SALMON=/home/rocesv/salmon-1.5.2_linux_x86_64/bin/
GROUP=$(ls ${INPUT})

for g in $GROUP
	do
	echo "Processing $g files ... | $(date)"
	LIB=$(ls ${INPUT}${g})
	for l in $LIB
		do
		echo "		Processing $g/$l files ... | $(date)"
		if [ "${l}" == "Paired" ]; then
		for f in ${INPUT}/${g}/${l}/*_1.fastq
			do
			base=${f%_1*}
			foward="${base}_1.fastq"
	  		reverse="${base}_2.fastq"
			echo " 			Processing $base files ... | $(date)"
			# script for paired
			mkdir ${OUTPUT}${g}/$(basename ${base})
			$SALMON/salmon quant --index $INDEX --libType ISR --dumpEq --hardFilter -1 ${foward} -2 ${reverse} -p 16 --output ${OUTPUT}${g}/$(basename ${base})/$(basename ${base}).out 
			done
		elif [ "${l}" == "Single" ]; then
		for f in ${INPUT}/${g}/${l}/*
			do
			echo " 			Processing  $f files ... | $(date)"
			# script for single
			mkdir ${OUTPUT}${g}/$(basename ${f})
			$SALMON/salmon quant --index $INDEX --libType A --dumpEq --hardFilter -r ${f} -p 16 --output ${OUTPUT}${g}/$(basename ${f})/$(basename ${f}).out
			done
		else
			# no script for unpaired
			echo "		CONTINUE"
			continue
		fi
		done
	done
