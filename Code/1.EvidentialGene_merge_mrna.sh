#! /bin/bash

# params

ncpu=16
maxmem=61000
evigene=/home/rocesv/evigene/scripts/
output=/mnt/g/Pra-GE-ATLAS/ATLAS/Pinus_radiata/1.Assemblies/EvidentialGene/1.output/Pinus_radiata_all/
input=/mnt/g/Pra-GE-ATLAS/ATLAS/Pinus_radiata/1.Assemblies/EvidentialGene/0.input/Pinus_radiata_all/Pinus_radiata_all_renamed.cdna

# checks

export PATH=/home/rocesv/cdhit/:$PATH
export fastanrdb=/home/rocesv/exonerate-2.2.0-x86_64/bin/fastanrdb
#blastn already in /usr/bin/ should not be a problem

# command

cd $output
echo $evigene/prot/tr2aacds.pl -tidy -NCPU $ncpu -MAXMEM $maxmem -log -cdna $input
$evigene/prot/tr2aacds.pl -tidy -NCPU $ncpu -MAXMEM $maxmem -log -cdna $input

# translation

transeq -sequence "${j}.allokay.cdna" -outseq "${j}.allokay.aa" -frame 6 -clean