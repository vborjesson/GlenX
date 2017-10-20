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

# make region specific bam file and convert to fasta file
samtools view -b $1 $2 $7 > $3_bam/$4.bam 
samtools view -h -F 2048 $3_bam/$4.bam | samtools view -ShF 1024 - | samtools view -ShF 512 - | samtools view -SF 256 - | awk '{print ">"$1"\n"$10}' > $3_fasta/$4.fasta
samtools view -h -F 2048 $3_bam/$4.bam | samtools view -ShF 1024 - | samtools view -ShF 512 - | samtools view -SF 256 - | grep "SA:Z:" > $3_fasta/$4_presplits.sam

# See if SVs in x_split.sam is in the correct region. 
python check_region.py --sam $3_fasta/$4_presplits.sam --region1 $2 --region2 $7 --out $3_fasta/$4_splits.sam
python consensus.py $3_fasta/$4_splits.sam $3_fasta > $3_fasta/consensus.fa

# make de novo assembly with three different kmer sizes; 30, 50 and 70
ABYSS -k 30 -c 3 -o $3_assembly/$4_30_contig.fa $3_fasta/$4.fasta
ABYSS -k 50 -c 3 -o $3_assembly/$4_50_contig.fa $3_fasta/$4.fasta
ABYSS -k 70 -c 3 -o $3_assembly/$4_70_contig.fa $3_fasta/$4.fasta
ABYSS -k 90 -c 3 -o $3_assembly/$4_90_contig.fa $3_fasta/$4.fasta
SSAKE -p 0 -w10 -f $3_fasta/$4.fasta -b $3_assembly/$4_sake
#mkdir $3_assembly/$4_velvet
#velveth $3_assembly/$4_velvet 31 -fasta -short $3_fasta/$4.fasta 
#velvetg $3_assembly/$4_velvet -cov_cutoff 5 

# rename id to be kmer specific in generated files and merge to one file
sed 's/^>/>30_/' $3_assembly/$4_30_contig.fa > $3_assembly/$4_merged.contig.fa
sed 's/^>/>50_/' $3_assembly/$4_50_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>70_/' $3_assembly/$4_70_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>90_/' $3_assembly/$4_90_contig.fa >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>sake_/' $3_assembly/$4_sake.contigs >> $3_assembly/$4_merged.contig.fa
sed 's/^>/>consensus_/' $3_fasta/consensus.fa >> $3_assembly/$4_merged.contig.fa
#sed 's/^>/>velvet_/' $3_assembly/$4_velvet/contigs.fa >> $3_assembly/$4_merged.contig.fa

# map back new contigs to reference genome
bwa mem -x intractg $6 $3_assembly/$4_merged.contig.fa > $3_assembly/$4_mapped.sam  


