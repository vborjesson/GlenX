#!/bin/bash -l

#SBATCH -A b2016296
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -J GlenX
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vanja.borjesson@gmail.com

module load bioinfo-tools
module load bwa
module load velvet
module load abyss
module load samtools

python GlenX.py --vcf ../GlenX_data/P2109_103_top10.vcf --bam ../GlenX_data/P2109_103.clean.dedup.recal.bam --tab ../GlenX_data/P2109_103.10000.tab --ID test01 --bwa_ref /proj/b2016296/private/nobackup/annotation/human_g1k_v37.fasta