#!bin/bash

INPUT=/mnt/g/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/00.Chunks/EviGene_allokay/prot/HomePC/
OUTPUT=/mnt/g/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/1.InterProScan5/EviGene_allokay/
IP5=/home/rocesv/my_interproscan/interproscan-5.52-86.0/interproscan.sh
TMP=/mnt/g/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/1.InterProScan5/0.tmp/

for f in $INPUT/*
	do
	base=$(basename $f)
	echo "Processing IP5 ${base} file ... | $(date)"
	bash $IP5 -cpu 18 -d $OUTPUT -f TSV,GFF3 -goterms -iprlookup -pa -t p -T $TMP -i $f
	done
