#!/bin/sh

# Make new region-specific bam-files, convert to fasta, Assembly and map contigs to genome generating a sam file. 

#$1	Path to bam-file  
#$2 Region of interest
#$3 bam file-name
#$4 fasta file-name 
#$5	sam file-name
#$6 fasta reference 
#$7 ID  
#$8 region 2

#module load bioinfo-tools
#module load samtools 

samtools view -b $1 $2 $8 > $3
samtools view $3 | awk '{print ">"$1"\n"$10}' > $4

bwa mem -x intractg 


