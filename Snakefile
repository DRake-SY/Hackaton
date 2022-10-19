samples=["SRR636533"]

rule all:
    input:
        expand("fastqc/{SAMPLE}_{n}_fastqc.html", SAMPLE=samples, n=[1,2])


rule fastqc:
    input:
        sample1="{sample}_1.fastq",
        sample2="{sample}_2.fastq"
    output:
        html1="fastqc/{sample}_1_fastqc.html",
        html2="fastqc/{sample}_2_fastqc.html"
    threads: 2

    shell: "fastqc -o fastqc -t {threads} {input.sample1} {input.sample2}"