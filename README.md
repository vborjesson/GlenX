# GlenX
Tool for better genotyping of whole genome sequencing (WGS) data. 

![alt text](https://github.com/vborjesson/GlenX/blob/master/Glen.png)

Takes an vcf, bam and tab-file with read depth per 100 bp as input. 

dependencies: 
	GlenX uses three different assembly tools for de novo assemly and mapping; abyss, velveth, SSAKE, BWA and samtools. All these softwares needs to be installed. 
	You will also need a SQlite-database with reference genome GC-content and mappability score / 100bp. Look at; Create database. 
```
ABySS
SSAKE
velvet
BWA
samtools
clustalw
GlenX_DB
consensus (https://github.com/J35P312/SplitVision) 
cdhit - (https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#CDHITEST)
```

## Create database
requirements; pandas, numpy
An database needs to be created in order for normalization of data before analyses can be done. Download mappability scores for reference genome hg19 from rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability and convert the bigwig file to readable bedGraph. Do also download the hg19 reference genome. 

```
python gc_map_array.py --hg19 /path/to/refgenome --bed /path/to/bedGraph
```
This will create a GlenX.db database in the GlenX directory

## Run GlenX
```
python GlenX --vcf --bam --tab --bwa_ref --norm_db
```
This will generate a vcf file with new predicted breakpoints, genotyping and SV classification. 

2017 - 10 - 23

GlenX_DB is done. 
Now im trying to figure out how to normalize the data from two databases. So far; I select what I want from the two databases where they are joined. 
According to this formula: 
	r_i_norm = r_i * m_all / m_gc 
	where:
	r_i_norm is the normaliized red coverage
	r_i is the cittent read coverage
	m_all is the average read_coverage
	m_gc is the read coverage for all regions with the same gc-content



