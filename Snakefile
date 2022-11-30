samples=["SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589"]


rule all:
    input:
        expand(["fastqc/{SAMPLE}_1_fastqc.html","fastqc/{SAMPLE}_2_fastqc.html","chromosome/ref.fa","chromosome/chr_annotation.gtf","star/{SAMPLE}.bam.bai","result.counts"], SAMPLE=samples)


rule download_fastq:
 output:
    "samples/{SAMPLE}_1.fastq","samples/{SAMPLE}_2.fastq"
 singularity:
  "docker://pegi3s/sratoolkit:2.10.0"
 threads: 16
 resources: load = 25 #a tester sans 
 shell:
  "fasterq-dump {wildcards.SAMPLE} --progress -O samples"


##Fastqc module
rule fastqc:
    input:
        sample1="samples/{SAMPLE}_1.fastq",
        sample2="samples/{SAMPLE}_2.fastq"
    output:
        html1="fastqc/{SAMPLE}_1_fastqc.html",
        html2="fastqc/{SAMPLE}_2_fastqc.html"
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
    resources: load=25
    shell: "wget -O chromosome/chr_annotation.gtf.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz\
    && gunzip chromosome/chr_annotation.gtf.gz"


#Index Star
rule index:
    input:"chromosome/ref.fa"
    output:"chromosome/SAindex"
    singularity:"docker://drakesy/hackaton:starv2"
    threads: 16
    resources: load=25
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
    resources: load=25
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
    resources: load=25
    shell: "samtools index {input}"

#feature_count
rule counting_reads:
 input:
  gtf="chromosome/chr_annotation.gtf"
 output:
  "result.counts"
 singularity:
  "docker://drakesy/hackaton:subread"
 threads: 4
 resources: load=25
 shell:
  """
  for samp in {samples};
  do 
   featureCounts -T {threads} -t gene -g gene_id -s 0 -a {input.gtf} -o {output} star/$samp.bam
  done
  """


