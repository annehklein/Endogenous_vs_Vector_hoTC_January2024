# Endogenous_vs_Vector_hoTC_January2024

# Endogenous vs Vector hotc Sam project

This project was realised for Samantha Ginn (sginn@cmri.org.au) from CMRI. Data comprise two PCR products for NGS which are comprised of two transcripts. The aim would be to determine the percentage of each transcript in the sample. One is endogenous and the other is introduced by NHEJ.

These are the two sequences that should be in the amplicon which is 311 bp. Paired ends PE250 overlap has been used fro sequencing.

Endogenous hOTC: 
ACCAAGCTGTTGCTGACAAATGATCCATTGGAAGCAGCGCATGGAGGCAATGTATTAATTACAGACACTTGGATAAGCATGGGACAAGAAGAGGAGAAGAAAAAGCGGCTCCAGGCTTTCCAAGGTTACCAGGTTACAATGAAGACTGCTAAAGTTGCTGCCTCTGACTGGACATTTTTACACTGCTTGCCCAGAAAGCCAGAAGAAGTGGATGATGAAGTCTTTTATTCTCCTCGATCACTAGTGTTCCCAGAGGCAGAAAACAGAAAGTGGACAATCATGGCTGTCATGGTGTCCCTGCTGACAGATTA

Vector encoded hOTC: 
ACCAAGCTGCTGCTGACAAATGACCCTCTGGAGGCAGCACACGGAGGAAACGTGCTGATCACCGATACATGGATCAGCATGGGCCAGGAGGAGGAGAAGAAGAAGCGGCTGCAGGCCTTCCAGGGCTACCAGGTGACCATGAAGACAGCCAAGGTGGCCGCCTCCGACTGGACCTTTCTGCACTGCCTGCCCAGAAAGCCTGAGGAGGTGGACGATGAGGTGTTCTACTCCCCCAGGTCTCTGGTGTTTCCTGAGGCCGAGAATCGCAAGTGGACCATCATGGCCGTGATGGTGTCCCTGCTGACAGATTA

Dataset consist in:
Pool-1-LFM11773_L1_1.fq Pool-1-LFM11773_L1_2.fq
Pool-2-LFM11774_L1_1.fq Pool-2-LFM11774_L1_2.fq


All data are stored [here](https://csiroau-my.sharepoint.com/:f:/r/personal/kle065_csiro_au/Documents/CMRI%20projects/Sam/January_2024?csf=1&web=1&e=KJTdst).

# Workflow

**Step 1: QC**

**FastQC output**
- [Pool1- Forward](https://github.com/annehklein/Endogenous_vs_Vector_hoTC_January2024/blob/main/Pool-1-LFM11773_L1_1_fastqc.html) / [Pool1- reverse](https://github.com/annehklein/Endogenous_vs_Vector_hoTC_January2024/blob/main/Pool-1-LFM11773_L1_2_fastqc.html)
- [Pool2- Forward](https://github.com/annehklein/Endogenous_vs_Vector_hoTC_January2024/blob/main/Pool-2-LFM11774_L1_1_fastqc.html)/ [Pool2 Reverse](https://github.com/annehklein/Endogenous_vs_Vector_hoTC_January2024/blob/main/Pool-2-LFM11774_L1_2_fastqc.html)


Quality filtering remove reads of quality under 20:

bbduk.sh in1=Pool1-LFE12852_L1_1.fq in2=Pool1-LFE12852_L1_1.fq out1=Pool1-LFE12852_L1_filt_1.fq out2=Pool1-LFE12852_L1_filt_2.fq maq=20
```
Pool 1
Input:                  	13402958 reads 		3350739500 bases.
Low quality discards:   	1006390 reads (7.51%) 	251597500 bases (7.51%)
Total Removed:          	1006390 reads (7.51%) 	251597500 bases (7.51%)
Result:                 	12396568 reads (92.49%) 	3099142000 bases (92.49%)

```
```
Pool 2
Input:                  	13662440 reads 		3415610000 bases.
Low quality discards:   	886138 reads (6.49%) 	221534500 bases (6.49%)
Total Removed:          	886138 reads (6.49%) 	221534500 bases (6.49%)
Result:                 	12776302 reads (93.51%) 	3194075500 bases (93.51%)
```

**Step 2: Merging forward and reverse reads**

```
Pool 1
Pairs:               	6388151
Joined:              	5295676   	82.898%
Ambiguous:           	298005   	4.665%
No Solution:         	794470   	12.437%
Too Short:           	0       	0.000%

Avg Insert:          	314.7
Standard Deviation:  	11.3
Mode:                	316

Insert range:        	94 - 491
90th percentile:     	316
75th percentile:     	316
50th percentile:     	316
25th percentile:     	316
10th percentile:     	313
```
```
Pool 2
Pairs:               	6198284
Joined:              	4045867   	65.274%
Ambiguous:           	705897   	11.389%
No Solution:         	1446520   	23.337%
Too Short:           	0       	0.000%

Avg Insert:          	315.6
Standard Deviation:  	8.1
Mode:                	316

Insert range:        	95 - 481
90th percentile:     	316
75th percentile:     	316
50th percentile:     	316
25th percentile:     	316
10th percentile:     	314
```


**Step 3: Split per barcode**

|**Pool 1**||||
|:-----------|:-----------|:-----------|:-----------|
|**Barcode**|**BC sequence** 	|**read count**|**percentage**|
|1|	GTTCA	|484907	|11.98|
|2	|GTCAT	|489370	|12.095|
|3|	CTGTA|	472124|	11.67|
|4	|GTATT|	401049|	9.91|
|5|	CTAGT|	421016	|10.41|
|6	|ACTTC|	416062	|13.58|
|7	|CCTAT|	444982|	10.998|
|8	|ACTGA|	493340|	12.19|
|9	|TCCAA|	293	|0.00724|
|10	|GCATT|	480|	0.01|
		
|**Pool 2**||||
|:-----------|:-----------|:-----------|:-----------|
|**Barcode**|**BC sequence** 	|**read count**|**percentage**|
|1|	GTTCA|	812588|	15.34|
|2	|GTCAT|	679785|	12.84|
|3|	CTGTA|	789781|	14.91|
|4|	GTATT|	671136|	12.67|
|5|	CTAGT|	721248	|13.62|
|6	|ACTTC	|79|	0.0015|
|7	|CCTAT|	93|	0.0017|
|8|	ACTGA|	148|	0.002|
|9	|TCCAA|	414465	|7.83|
|10	|GCATT	|541358|	10.22|

				


**Step 4: Align reads against Endogenous hOTC and Vector encoded hOTC**

Alignement using bwa.

The two reference sequences are of similar lentgh and have an identity of 19% so all reads align agaisnt both reference.

Need to add an extra step we only keep alignemnt with less than 5 different nucleotides

samtools view -H {input.r1} > {params.temp_header}
samtools view {input.r1} | awk -F '\t' '{{for(i=1;i<=NF;i++){{if($i ~ /^NM:i:/){{split($i,a,":"); if(a[3] <= 5){{print}}}}}}}}' >> {params.temp_header}
samtools view -bS {params.temp_header} > {output.filtered_bam}

|**Pool 1**|||||||
|:-----------|:-----------|:-----------|:-----------|:-----------|:-----------|:-----------|
|**Barcode**|**BC sequence** 	|**Endogenous reads**	|	**% Endogenous reads**	|**Vector reads**	|	**% Vector reads**	|
|1	|GTTCA|			479576	|98.9		|327	|0.067|
|2|	GTCAT|			485117|	99.13	|	328|	0.067|
|3|	CTGTA|			466857|	98.88|		293|	0.06|
|4	|GTATT|			390549	|97.38	|	262	|0.06|
|5	|CTAGT|		407200|	96.72		|676|	0.16|
|6	|ACTTC|			412262|	99.09	|	123|	0.029|
|7	|CCTAT|	432029|	97.089|		369|	0.08|
|8|	ACTGA|	453105|	91.84	|	145|	0.029|
|9	|TCCAA|	272|92.83	|	8|	2.73|
|10	|GCATT	|		462|	96.25	|	0|	0|

								
|**Pool 2**||||||
|:-----------|:-----------|:-----------|:-----------|:-----------|:-----------|
|**Barcode**|**BC sequence** | **Endogenous reads**	|	**% Endogenous reads**	|**Vector reads**	|	**% Vector reads**	|
|1|	GTTCA|	191945|23.62	|	614769|	75.65|
|2|	GTCAT	|139113|	20.46	|526180|	77.4|
|3|	CTGTA|	218470	|27.66	|	81089|	10.26|
|4|	GTATT|	308790|46	|	352743|	52.55|
|5	|CTAGT|		250269|	34.69		|453264|	62.84|
|6	|ACTTC|			70|	88.60	|	4|	5.06|
|7	|CCTAT|	86|	92.47|		1|1.075|
|8|	ACTGA|	125	|84.45	|6|	4.05|
|9	|TCCAA|	381796|	92.12	|	14838|	3.58|
|10|	GCATT	|526933|97.33	|1414	|0.26|
										


