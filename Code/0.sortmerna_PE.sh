#!bin/bash
FILES=/mnt/g/ATLAS_Pra/ATLAS/fastq_0/PE/*_1.fastq

for f in $FILES
do
 # log FILE-date  
 base=${f%_1*}
 foward="${base}_1.fastq"
 reverse="${base}_2.fastq"
 echo "Processing $base files ... | $(date)"
 
 # PE sortmerna execution: try it alone:other names format
 sortmerna -ref /mnt/g/ATLAS_Pra/tools/sortmerna/db/rfam-5s-database-id98.fasta -ref /mnt/g/ATLAS_Pra/tools/sortmerna/db/silva-arc-23s-id98.fasta -ref /mnt/g/ATLAS_Pra/tools/sortmerna/db/silva-bac-23s-id98.fasta -ref /mnt/g/ATLAS_Pra/tools/sortmerna/db/silva-euk-28s-id98.fasta -ref /mnt/g/ATLAS_Pra/tools/sortmerna/db/rfam-5.8s-database-id98.fasta -ref /mnt/g/ATLAS_Pra/tools/sortmerna/db/silva-arc-16s-id95.fasta -ref /mnt/g/ATLAS_Pra/tools/sortmerna/db/silva-bac-16s-id90.fasta -ref /mnt/g/ATLAS_Pra/tools/sortmerna/db/silva-euk-18s-id95.fasta -reads $foward -reads $reverse --fastx --workdir /mnt/g/ATLAS_Pra/tools/sortmerna/run/ --other --paired_in True -v True --threads 5 --num_alignments 1 -m 57000 -out2
 # Move files to desired directory & rename
 mv /mnt/g/ATLAS_Pra/tools/sortmerna/run/out/other_fwd.fastq /mnt/g/ATLAS_Pra/ATLAS/fastq_sortmerna/PE/$(basename $foward)
 mv /mnt/g/ATLAS_Pra/tools/sortmerna/run/out/other_rev.fastq /mnt/g/ATLAS_Pra/ATLAS/fastq_sortmerna/PE/$(basename $reverse)
 # Remove files and prepare next run
 rm -rf /mnt/g/ATLAS_Pra/tools/sortmerna/run/kvdb/*
 rm -rf /mnt/g/ATLAS_Pra/tools/sortmerna/run/out/*
done


