samples="SRR636533"


rule all:
    input:
        expand("fastqc/{SAMPLE}_{n}_fastqc.html", SAMPLE=samples, n=[1,2])


#can't download the sra file with prefetch

##fastq-dump
rule download_fastq:
    input:
        sra = "{samples}.sra"
    output: "samples/{samples}_{n}.fastq"
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
        "docker://Docker/fastqc/hackaton:fastqc"
        #Don't know if the code take it

    shell: "fastqc -o fastqc -t {threads} {input.sample1} {input.sample2}"


##Download chromosome


