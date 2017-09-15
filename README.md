# GlenX
Tool for better genotyping of whole genome sequencing (WGS) data. 

![alt text](https://github.com/vborjesson/GlenX/blob/master/Glen.png)

### Startup project: 
Try to find a good de novo assembly tool. I try Velvet, abyss and spades. To align the assembled reads (scaffolds) back to the ref genome again, I use bwa mem. 

#### What I tried so far: 

Get reads from specific region (start +/- 1kb - end +/- 1kb): 

```
samtools view -b file.bam "chrX:XXXXXX-XXXXX" "chrX:XXXXXX-XXXXXX" > regionX.bam
samtools view regionX.bam | awk '{print ">"$1"\n"$10}' > regionX.fasta

velveth velveth_out regionX.fasta
abyss-pe k-50 name=name in=regionX.fasta
spades.py --12 regionX.bam -o spades_out 

bwa mem -x intractg ref.fa assembled.fasta > aligned_assambled.sam
```
Velveth returns a fasta-file with all reads sorted, and after bwa I get a long sam-file - multiple steps -> i think I will skip this one!
Abyss - returns nice contigs. One good Ide could be to map contigs from different kmers. But; how do I get statistics from this? Look at sam! 
spades - finns på UPPMAX, ska testköra med inv_
SSAKE - Returnerar fler och kortare contigs -> ska kolla igenom!
tasr - ? 

Map reads back to contigs to get stats. - skip! This will not generate the maps that was originally created. 

Try to find other ways to get statistics for support!

Genotyping and classify type of SV;
From the mapped contigs a sam-file is generated. sam-file gives us information about direction of the reads and if they are matched or not matched to reference genome. 
Looking at the NOT matched contigs (S), I will get information about a SV and where it starts and ends. If this is mapped inside our regions of interest, we will call it a SV. And if we also have a contig that maps as match (M) to the reference genome, we can classify it as a heterozygous SV.   

Thoughts; 
go back to fasta-file generated by assembly tool, information about statistics? 

2017 - 09 -13
	A few changes most be done on GlenX. 
		1.	Read in TIDDIT tab-file with read coverage for every 100th bp (default is 500, is that enough?).
		2.	Count average coverage (exclude x- and y-chromosome due to different coverage?)
		3. 	Update counting breakpoints adding D as important information. Take into count reverse strand (flip cigar). 
		4. 	Make new array saving soft clips; strand, chr, pos, score, strand, chr, pos, score, cigar, contig length, number of alternative pos)
		5. 	Save match as before in array!
		6. 	Compare the two arrays, if not in region discard! 
		7. 	If more than 1 in array; make some kind of scoring, selecting the most important one. For ex. a lot of cigars (complex), very short contig, 	a lot of alternative positions, low score etc.. Best score stays!
		8. 	Classify using tab-file, strand, position etc. 
		9. 	Write to new vcf-file! 

	
	1.

	After last meeting with J, we will run ABYSS like this: 

	ABYSS -k 30 -c X -o x_30 fasta
	ABYSS -k 50 -c X -0 x_50 fasta  
	ABYSS -k 70 -c X -o x_70 fasta
	Where X is the read coverage divided by 4.

	cat x_* >> merged_inv.fa

	bwa mem -x intractg /proj/b2016296/private/nobackup/annotation/human_g1k_v37.fasta merged_inv.fa > aligned_merged_inv.sam

	Run GlenX!

	2.

	SSAKE -f fasta -p0 -w1 -b sake_inv
	bwa mem -x intractg /proj/b2016296/private/nobackup/annotation/human_g1k_v37.fasta sake_inv.contig > aligned_sake_inv.sam

	SSAKE generates very short but many contigs. Can not alternative alignment (end_positions).
	After discussion with Jesper, we will skip SSAKE. 

	3. 
	Have already tried with abyss k-55 this results in 
	chrA 11 pos 30895551 chrB 11 pos 26504326 INV 0/1 -> pretty close to the expected chrA 11 pos 30895624 chrB 11 pos 26504276 INV 0/1

	Updated GlenX taking hard clipping into count, this will change the breakpoints. 


  


