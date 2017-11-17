#!/bin/sh

# Make new region-specific bam-files, convert to fasta, Assembly and map contigs to genome generating a sam file. 

#$1	Path to bam-file  
#$2 Region of interest
#$3 ID
#$4 region ID 
#$5 minimum coverage accepted 
#$6 bwa reference genome
#$7 region 2

#module load bioinfo-tools
#module load samtools 
#module load abyss
#module load bwa
#module load clustalw

# create a new folder for consensus output 
mkdir $3_consensus

# make region specific bam file and convert to fasta file
samtools view -b $1 $2 $7 > $3_bam/$4.bam 
samtools view -h -F 2048 $3_bam/$4.bam | samtools view -ShF 1024 - | samtools view -ShF 512 - | samtools view -SF 256 - | awk '{print ">"$1"\n"$10}' > $3_fasta/$4.fasta
samtools view -h -F 2048 $3_bam/$4.bam | samtools view -ShF 1024 - | samtools view -ShF 512 - | samtools view -SF 256 - | grep "SA:Z:" > $3_consensus/$4_presplits.sam

# See if SVs in x_presplit.sam is in the correct region. This is a prefilter before consensus.
python modules/check_region.py --sam $3_consensus/$4_presplits.sam --region1 $2 --region2 $7 --out $3_consensus/$4_splits.sam
mkdir $3_consensus/$4_consensus
python consensus.py $3_consensus/$4_consensus $3_consensus/$4_splits.sam > $3_consensus/$4_consensus.fa

# make de novo assembly with three different kmer sizes; 30, 50 and 70
ABYSS -k 30 -c 3 -o $3_assembly/$4_30_1_contig.fa $3_fasta/$4.fasta
ABYSS -k 30 -t 10 -c 3 -o $3_assembly/$4_30_2_contig.fa $3_fasta/$4.fasta
ABYSS -k 30 -t 30 -c 3 -o $3_assembly/$4_30_3_contig.fa $3_fasta/$4.fasta
ABYSS -k 50 -c 3 -o $3_assembly/$4_50_1_contig.fa $3_fasta/$4.fasta
ABYSS -k 50 -t 10 -c 3 -o $3_assembly/$4_50_2_contig.fa $3_fasta/$4.fasta
ABYSS -k 50 -t 30 -c 3 -o $3_assembly/$4_50_3_contig.fa $3_fasta/$4.fasta
ABYSS -k 70 -c 3 -o $3_assembly/$4_70_1_contig.fa $3_fasta/$4.fasta
ABYSS -k 70 -t 10 -c 3 -o $3_assembly/$4_70_2_contig.fa $3_fasta/$4.fasta
ABYSS -k 70 -t 30 -c 3 -o $3_assembly/$4_70_3_contig.fa $3_fasta/$4.fasta
ABYSS -k 90 -c 3 -o $3_assembly/$4_90_1_contig.fa $3_fasta/$4.fasta
ABYSS -k 90 -t 10 -c 3 -o $3_assembly/$4_90_2_contig.fa $3_fasta/$4.fasta
ABYSS -k 90 -t 30 -c 3 -o $3_assembly/$4_90_3_contig.fa $3_fasta/$4.fasta
timeout 1m SSAKE -p 0 -w10 -f $3_fasta/$4.fasta -b $3_assembly/$4_sake
mkdir $3_assembly/$4_velvet
velveth $3_assembly/$4_velvet 31 -fasta -short $3_fasta/$4.fasta 
velvetg $3_assembly/$4_velvet -cov_cutoff 5 

# rename id to be kmer specific in generated files and merge to one files
sed 's/^>/>30_1_/' $3_assembly/$4_30_1_contig.fa > $3_assembly/$4_merged.contig.fa
sed 's/^>/>30_2_/' $3_assembly/$4_30_2_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>30_3_/' $3_assembly/$4_30_3_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>50_1_/' $3_assembly/$4_50_1_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>50_2_/' $3_assembly/$4_50_2_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>50_3_/' $3_assembly/$4_50_3_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>70_1_/' $3_assembly/$4_70_1_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>70_2_/' $3_assembly/$4_70_2_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>70_3_/' $3_assembly/$4_70_3_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>90_1_/' $3_assembly/$4_90_1_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>90_2_/' $3_assembly/$4_90_2_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>90_3_/' $3_assembly/$4_90_3_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>sake_/' $3_assembly/$4_sake.contigs >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>consensus_/' $3_consensus/$4_consensus.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>velvet_/' $3_assembly/$4_velvet/contigs.fa >> $3_assembly/$4_merged.contig.fa

# map back new contigs to reference genome
bwa mem -x intractg $6 $3_assembly/$4_merged.contig.fa > $3_assembly/$4_mapped.sam  


