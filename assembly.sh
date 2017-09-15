#!/bin/sh

# Make new region-specific bam-files, convert to fasta, Assembly and map contigs to genome generating a sam file. 

#$1	Path to bam-file  
#$2 Region of interest
#$3 bam file-name
#$4 fasta file-name 
#$5	sam file-name
#$6 assembly map
#$7 fasta reference 
#$8 ID  
#$9 region 2

#module load bioinfo-tools
#module load samtools 

samtools view -b $1 $2 $9 > $3
samtools view $3 | awk '{print ">"$1"\n"$10}' > $4
abyss-pe k-55 name= $8 in= $4 
bwa mem -x intractg 


