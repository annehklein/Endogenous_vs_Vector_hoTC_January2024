barcodes = ["GTTCA","GTCAT","CTGTA","GTATT","CTAGT","ACTTC","CCTAT","ACTGA","TCCAA","GCATT"]
read_qual=20
rule all:
    input:expand('reads/Pool{pool}_L1_1.fq', pool={"1-LFE12852","2-LFE12853","3-LFE12854","4-LFE12855"}),
          expand('reads/Pool{pool}_L1_2.fq', pool={"1-LFE12852","2-LFE12853","3-LFE12854","4-LFE12855"}),
          #expand('data/Pool{pool}_{number}.filt.fastq',pool={"1-LFE12852","2-LFE12853","3-LFE12854","4-LFE12855"},number={"1","2"}),
          expand('data/Pool{sample}_1.filt.fastq',sample={"1-LFE12852","2-LFE12853","3-LFE12854","4-LFE12855"}),
	  expand('data/Pool{sample}_2.filt.fastq',sample={"1-LFE12852","2-LFE12853","3-LFE12854","4-LFE12855"}),        	       expand('data/Pool{sample}.fq',sample={"1-LFE12852","2-LFE12853","3-LFE12854","4-LFE12855"}),
	  expand('data/Pool{sample}.{seq}.fq',sample={"1-LFE12852","2-LFE12853","3-LFE12854","4-LFE12855"},seq=barcodes),
	  expand('data/Pool{sample}.{seq}.log',sample={"1-LFE12852","2-LFE12853","3-LFE12854","4-LFE12855"},seq=barcodes),
          expand('data/Endogenous{sample}.{seq}.bam',sample={"1-LFE12852","2-LFE12853","3-LFE12854","4-LFE12855"},seq=barcodes),
	  expand('data/Vector{sample}.{seq}.bam',sample={"1-LFE12852","2-LFE12853","3-LFE12854","4-LFE12855"},seq=barcodes),
	  expand('data/{reference}{sample}.{seq}.bam',reference={"Endogenous","Vector"},sample={"1-LFE12852","2-LFE12853","3-LFE12854","4-LFE12855"},seq=barcodes),
	  expand('data/Filtered_{reference}{sample}.{seq}.bam',reference={"Endogenous","Vector"},sample={"1-LFE12852","2-LFE12853","3-LFE12854","4-LFE12855"},seq=barcodes),
#	  expand('data/header.{sample}.sam',sample={"1-LFE12852","2-LFE12853","3-LFE12854","4-LFE12855"})
        #  expand('data/S1-LEM4848_L1_{nombre}.filt.{seq}.bam.bai',nombre={"1","2"}, seq=barcodes),
         # expand('data/S1-LEM4848_L1_{nombre}.{seq}.{allele}.bam',nombre={"1","2"}, seq=barcodes,allele={"ref.G","alt.C"}),
        #  expand('data/S1-LEM4848_L1_{nombre}.{seq}.{allele}.bam.bai',nombre={"1","2"}, seq=barcodes,allele={"ref.G","alt.C"})
          #expand('output_OTC5_A_C/S1-LEM4848_L1_{nombre}.{seq}.{allele}.fq',nombre={"1","2"}, seq=barcodes,allele={"ref.C","alt.A"})
          #expand('reads/OTC3-LEB3139_L1_{nombre}.{seq}.{allele}.fq',nombre={"1","2"}, seq=barcodes,allele={"ref.A","alt.G"})


rule qual_filt:
    input:
        r1="reads/Pool{pool}_L1_1.fq",
        r2="reads/Pool{pool}_L1_2.fq"
    output:
        out1="data/Pool{pool}_1.filt.fastq",
        out2="data/Pool{pool}_2.filt.fastq"
    threads: 5
    shell:
        "bbduk.sh in1={input.r1} in2={input.r2} out1={output.out1} out2={output.out2} maq={read_qual}"

rule merge:
	input:
		r1='data/Pool{sample}_1.filt.fastq',
		r2='data/Pool{sample}_2.filt.fastq' 
	output:
		merged='data/Pool{sample}.fq'

	shell:
		"""
		bbmerge.sh in1={input.r1} in2={input.r2} \
		out={output.merged}
		"""
rule barcode_filt:
    input:
        "data/Pool{sample}.fq"
    output:
       filt = "data/Pool{sample}.{seq}.fq",
       log = "data/Pool{sample}.{seq}.log"
    threads: 5
    params:
       seq= "{seq}"
    shell:
        """
        python3 filter.py -i {input} -o {output.filt} -seq {params.seq} -log {output.log}
        """

rule align_endogenous:
    input:
        r1 = "data/Pool{sample}.{seq}.fq",
        index = expand('reference/Endogenous_hOTC{ext}', ext=['.amb', '.ann', '.bwt', '.pac', '.sa'])
    output:
        bam = 'data/Endogenous{sample}.{seq}.bam'
    params:
        prefix = lambda wildcards, input: ".".join(input.index[0].split('.')[:-1])
    shell:
        """
        bwa mem -M -t 2 -p  {params.prefix} {input.r1}   2> bwa.err | samtools sort -o {output.bam} 
        """
rule align_vector:
    input:
        r1 = "data/Pool{sample}.{seq}.fq",
        index = expand('reference/Vector_encoded_hOTC{ext}', ext=['.amb', '.ann', '.bwt', '.pac', '.sa'])
    output:
        bam = 'data/Vector{sample}.{seq}.bam'
    params:
        prefix = lambda wildcards, input: ".".join(input.index[0].split('.')[:-1])
    shell:
        """
        bwa mem -M -t 2 -p  {params.prefix} {input.r1}   2> bwa.err | samtools sort -o {output.bam} 
        """
rule filter_aligment:
    input:
        r1 = "data/{reference}{sample}.{seq}.bam",
    output:
        filtered_bam = "data/Filtered_{reference}{sample}.{seq}.bam"
    params:
        temp_header = "data/header.{sample}.sam"
    shell:
        """
	samtools view -H {input.r1} > {params.temp_header}
        samtools view {input.r1} | awk -F '\t' '{{for(i=1;i<=NF;i++){{if($i ~ /^NM:i:/){{split($i,a,":"); if(a[3] <= 5){{print}}}}}}}}' >> {params.temp_header}
        samtools view -bS {params.temp_header} > {output.filtered_bam}
        rm {params.temp_header}
	"""
