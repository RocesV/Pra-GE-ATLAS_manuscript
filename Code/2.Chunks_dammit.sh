#!bin/bash

INPUT=/mnt/g/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/00.Chunks/Lace_allokay/nt/
OUTPUT=/mnt/g/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/2.dammit/Lace_okayall/
EPIOMIDES=/mnt/e/ATLAS_Pra/ATLAS/DB/EpiomidesGymnosperms.fa
DB=/mnt/e/ATLAS_Pra/ATLAS/DB/databases_/

for f in $INPUT/*
	do
	base=$(basename $f)
	echo "Processing dammit ${base} file ... | $(date)"
	mkdir $OUTPUT/$base
	dammit annotate $f -o $OUTPUT/$base/ --database-dir $DB --busco-group embryophyta --n_threads 18 --full --no-rename
	done
