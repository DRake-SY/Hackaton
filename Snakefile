samples=["SRR636533"]
chromosome=["2"]


rule all:
    input:
        expand(["fastqc/{SAMPLE}_{n}_fastqc.html","chromosome/ch{CHROMO}.fa"], SAMPLE=samples,CHROMO=chromosome, n=[1,2])


rule prefetch:
    output:"samples/{SAMPLE}.sra"
    shell: 'prefetch {wildcards.SAMPLE} -O samples \
    && mv samples/{wildcards.SAMPLE}/{wildcards.SAMPLE}.sra samples \
    && rm -r samples/{wildcards.SAMPLE}'

##fastq-dump
rule fastq_dump:
    input: 
        sra="samples/{SAMPLE}.sra"
    output: "samples/{SAMPLE}_{n}.fastq"
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

    #containerized:
        #"Docker/fastqc/dockerfile"
        #Don't know if the code take it

    shell: "fastqc -o fastqc -t {threads} {input.sample1} {input.sample2}"


##Download chromosome index


rule index:
    output: "chromosome/ch{CHROMO}.fa"

    shell: "wget -O chromosome/ch{wildcards.CHROMO}.fa.gz  https://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{wildcards.CHROMO}.fa.gz \
    && gunzip chromosome/ch{wildcards.CHROMO}.fa.gz \
    && STAR --runMode genomeGenerate \
    --genomeDir  chromosome/ \
    --genomeFastaFiles chromosome/ch{wildcards.CHROMO}.fa "


