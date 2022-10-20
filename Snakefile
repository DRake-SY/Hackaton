samples=["SRR636533"]

rule all:
    input:
        expand("fastqc/{SAMPLE}_{n}_fastqc.html", SAMPLE=samples, n=[1,2])


##Download fastq data



##Fastqc module
rule fastqc:
    input:
        sample1="{sample}_1.fastq",
        sample2="{sample}_2.fastq"
    output:
        html1="fastqc/{sample}_1_fastqc.html",
        html2="fastqc/{sample}_2_fastqc.html"
    threads: 2

    singularity:
        "docker://Docker/fastqc/hackaton:fastqc"
        #Don't know if the code take it

    shell: "fastqc -o fastqc -t {threads} {input.sample1} {input.sample2}"


##Download chromosome

