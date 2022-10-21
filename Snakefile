samples=["SRR636533"]


rule all:
    input:
        expand("fastqc/{SAMPLE}_{n}_fastqc.html", SAMPLE=samples, n=[1,2])




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

    container:
        "docker://Docker/fastqc/hackaton:fastqc"
        #Don't know if the code take it

    shell: "fastqc -o fastqc -t {threads} {input.sample1} {input.sample2}"


##Download chromosome
