#!bin/bash
FILES=/mnt/e/ATLAS_Pra/ATLAS/fastq_0/SE/*

for f in $FILES
do 
 # log FILE-date
 echo "Processing $(basename $f) file... | $(date)"
 
 # SE sortmerna execution
 sortmerna -ref /mnt/e/ATLAS_Pra/tools/sortmerna/db/rfam-5s-database-id98.fasta -ref /mnt/e/ATLAS_Pra/tools/sortmerna/db/silva-arc-23s-id98.fasta -ref /mnt/e/ATLAS_Pra/tools/sortmerna/db/silva-bac-23s-id98.fasta -ref /mnt/e/ATLAS_Pra/tools/sortmerna/db/silva-euk-28s-id98.fasta -ref /mnt/e/ATLAS_Pra/tools/sortmerna/db/rfam-5.8s-database-id98.fasta -ref /mnt/e/ATLAS_Pra/tools/sortmerna/db/silva-arc-16s-id95.fasta -ref /mnt/e/ATLAS_Pra/tools/sortmerna/db/silva-bac-16s-id90.fasta -ref /mnt/e/ATLAS_Pra/tools/sortmerna/db/silva-euk-18s-id95.fasta -reads $f --fastx --workdir /mnt/e/ATLAS_Pra/tools/sortmerna/run/ --other -v True --threads 5 --num_alignments 1 -m 29000
 # Move files to desired directory & rename
 mv /mnt/e/ATLAS_Pra/tools/sortmerna/run/out/other.fastq /mnt/e/ATLAS_Pra/ATLAS/fastq_sortmerna/SE/$(basename $f)
 # Remove files and prepare next run
 rm -rf /mnt/e/ATLAS_Pra/tools/sortmerna/run/kvdb/*
 rm -rf /mnt/e/ATLAS_Pra/tools/sortmerna/run/out/* 
done
