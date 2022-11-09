samples=["SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589"]
strands = ['0', '1', '2']

rule all:
    input:
        expand(["fastqc/{SAMPLE}_{n}_fastqc.html","chromosome/ref.fa", "chromosome/chr_annotation.gtf", "chromosome/SAindex","star/{SAMPLE}.bam.bai","gene_{SAMPLE}_strand_{STRAND}.counts"], SAMPLE=samples, STRAND=strands, n=[1,2])

rule prefetch:
    output:"samples/{SAMPLE}.sra"
    singularity: "docker://evolbioinfo/sratoolkit:v2.10.8"
    shell: 'vdb-config -i \
    && prefetch {wildcards.SAMPLE} -O samples \
    && mv samples/{wildcards.SAMPLE}/{wildcards.SAMPLE}.sra samples \
    && rm -r samples/{wildcards.SAMPLE}'

##fastq-dump
rule fastq_dump:
    input: 
        sra="samples/{SAMPLE}.sra"
    output: "samples/{SAMPLE}_{n}.fastq"

    singularity: "docker://biocontainers/sra-toolkit:v2.9.3dfsg-1b1-deb_cv1"
    shell: 'fastq-dump --split-files {input.sra} -O samples'
    

##Fastqc module
rule fastqc:
    input:
        sample1="samples/{sample}_1.fastq",
        sample2="samples/{sample}_2.fastq"
    output:
        html1="fastqc/{sample}_1_fastqc.html",
        html2="fastqc/{sample}_2_fastqc.html"
    threads: 2
    singularity:
        "docker://drakesy/hackaton:fastqcv2"
    shell: "fastqc -o fastqc -t {threads} {input.sample1} {input.sample2}"


##Download chromosome index
rule download_chromosome:
   output:
    "chromosome/ref.fa"
   shell:
    """
    for chromosome in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "MT" "X" "Y";
    do
        wget -O "$chromosome".fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome."$chromosome".fa.gz
    done
    gunzip -c *.fa.gz > {output} \
    && rm *.fa.gz
    """

#Genome annotation
rule download_genome_annotation:
    output: "chromosome/chr_annotation.gtf"

    shell: "wget -O chromosome/chr_annotation.gtf.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz\
    && gunzip chromosome/chr_annotation.gtf.gz"


#Index Star
rule index:
    input:"chromosome/ref.fa"
    output:"chromosome/SAindex"
    singularity:"docker://drakesy/hackaton:starv2"
    threads: 16
    shell:
     "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir chromosome/ --genomeFastaFiles {input}"


#mapping
rule mapping:
    input:
        "chromosome/SAindex",
        sample1="samples/{SAMPLE}_1.fastq",
        sample2="samples/{SAMPLE}_2.fastq"
    output:"star/{SAMPLE}.bam"
    singularity: "docker://drakesy/hackaton:starv2"
    threads: 16
    shell:"STAR --outSAMstrandField intronMotif \
--outFilterMismatchNmax 4 \
--outFilterMultimapNmax 10 \
--genomeDir chromosome \
--readFilesIn {input.sample1} {input.sample2} \
--runThreadN {threads} \
--outSAMunmapped None \
--outSAMtype BAM SortedByCoordinate \
--outStd BAM_SortedByCoordinate \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM 31000000000 \
> {output}"

#samtools
rule samtools:
    input: "star/{SAMPLE}.bam"
    output: "star/{SAMPLE}.bam.bai"
    singularity: "docker://drakesy/hackaton:samtools"
    shell: "samtools index {input}"


#couting_read
rule couting_reads:
    input:
        gtf="chromosome/chr_annotation.gtf",
        bam="star/{SAMPLE}.bam",
        bai="star/{SAMPLE}.bam.bai"
    output:
        "gene_{SAMPLE}_strand_{STRAND}.counts"
    singularity: "docker://drakesy/hackaton:subread"
    threads: 4
    shell:
     "featureCounts -T {threads} -t gene -g gene_id -s {wildcards.STRAND} -a {input.gtf} -o {output} {input.bam}"

