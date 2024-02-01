#!/usr/bin/env bash
# The name to show in queue lists for this job:
#SBATCH -J rnaspades_megagametophyte_2

# Number of desired cpus:
#SBATCH --ntasks=20

# Amount of RAM needed for this job:
#SBATCH --mem=1024gb

# The time the job will be running:
#SBATCH --time=7-00:00:00

# If you need nodes with special features uncomment the desired constraint line:
##SBATCH --constraint=slim
#SBATCH --constraint=bigmem
##SBATCH --constraint=cal

# Set output and error files
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

# To load some software (you can show the list with 'module avail'):
module load spades/3.14

# the program to execute with its parameters:
time /mnt/home/soft/spades/programs/x86_64/SPAdes-3.14.0-Linux/bin/spades.py --rna --checkpoints last --pe1-1 /mnt/home/users/bio_001_uovi/lvalledor/Paired_SRR13823451.sra_1.cor.fq.gz --pe1-2 /mnt/home/users/bio_001_uovi/lvalledor/Paired_SRR13823451.sra_2.cor.fq.gz -m 1000 --pe1-rf -o /mnt/scratch/users/bio_001_uovi/lvalledor/rnaspades_megagametophyte_2/
