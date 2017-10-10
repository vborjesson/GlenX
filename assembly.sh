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

# make region specific bam file and convert to fasta file
samtools view -b $1 $2 $9 > $3_bam/$4.bam 
samtools view $3_bam/$4.bam | awk '{print ">"$1"\n"$10}' > $3_fasta/$4.fasta

# make de novo assembly with three different kmer sizes; 30, 50 and 70
ABYSS -k 30 -c $5 -o $3_assembly/$4_30_contig.fa $3_fasta/$4.fasta
ABYSS -k 50 -c $5 -o $3_assembly/$4_50_contig.fa $3_fasta/$4.fasta
ABYSS -k 70 -c $5 -o $3_assembly/$4_70_contig.fa $3_fasta/$4.fasta
SSAKE -p 0 -w1 -f $3_fasta/$4.fasta -b $3_assembly/$4_sake
#mkdir $3_assembly/$4_velvet
#velveth $3_assembly/$4_velvet 31 -fasta -short $3_fasta/$4.fasta 
#velvetg $3_assembly/$4_velvet -cov_cutoff 5 

# rename id to be kmer specific in generated files and merge to one file
sed 's/^>/>30_/' $3_assembly/$4_30_contig.fa > $3_assembly/$4_merged.contig.fa
sed 's/^>/>50_/' $3_assembly/$4_50_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>70_/' $3_assembly/$4_70_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>sake_/' $3_assembly/$4_sake.contigs >> $3_assembly/$4_merged.contig.fa
#sed 's/^>/>velvet_/' $3_assembly/$4_velvet/contigs.fa >> $3_assembly/$4_merged.contig.fa

# map back new contigs to reference genome
bwa mem -x intractg $6 $3_assembly/$4_merged.contig.fa > $3_assembly/$4_mapped.sam  


