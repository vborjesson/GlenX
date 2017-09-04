# GlenX
Tool for better genotyping of whole genome sequencing (WGS) data. 

![alt text](https://github.com/vborjesson/GlenX/blob/master/Glen.png)

### Startup project: 
Try to find a good de novo assembly tool. I try Velvet, abyss and spades. To align the assembled reads (scaffolds) back to the ref genome again, I use bwa mem. 

#### What I tried so far: 

Get reads from specific region (start +/- 1kb - end +/- 1kb): 

```
samtools view -b "chrX:XXXXXX-XXXXX" "chrX:XXXXXX-XXXXXX" > regionX.bam
samtools view regionX.bam | awk '{print ">"$1"\n"$10}' > regionX.fasta

velveth velveth_out regionX.fasta
abyss-pe k-50 name=name in=regionX.fasta
spades.py --12 regionX.bam -o spades_out 

bwa mem -x intractg ref.fa assembled.fasta > aligned_assambled.sam
```
Velveth returns a fasta-file with all reads sorted, and after bwa I get a long sam-file - multiple steps -> i think I will skipp this one!
Abyss - returns nice contigs. One good Ide could be to map contigs from different kmers. But; how do I get statistics from this? Look at sam! 
spades - finns på UPPMAX, ska testköra med inv_
SSAKE - Returnerar fler och kortare contigs -> ska kolla igenom!
tasr - ? 

Map reads back to contigs to get stats. - skipp! This will not generate the maps that was originally created. 

Try to find other ways to get statistics!



