#!/usr/bin/env bash
# The name to show in queue lists for this job:
#SBATCH -J trinity_megagametophyte

# Number of desired cpus:
#SBATCH --ntasks=20

# Amount of RAM needed for this job:
#SBATCH --mem=128gb

# The time the job will be running:
#SBATCH --time=7-00:00:00


# If you need nodes with special features uncomment the desired constraint line:
##SBATCH --constraint=slim
##SBATCH --constraint=bigmem
##SBATCH --constraint=cal

# Set output and error files
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

# To load some software (you can show the list with 'module avail'):
module load trinity/2.10.0

# the program to execute with its parameters:
time trinity --seqType fq --max_memory 125G --CPU 16 --min_contig_length 350 --no_normalize_reads --bflyCPU 12 --left /mnt/home/users/bio_001_uovi/lvalledor/Paired_SRR13823451.sra_1.cor.fq.gz --right /mnt/home/users/bio_001_uovi/lvalledor/Paired_SRR13823451.sra_2.cor.fq.gz --output /mnt/scratch/users/bio_001_uovi/lvalledor/trinity_out_dir_megagametophyte/ --NO_SEQTK --SS_lib_type RF 
