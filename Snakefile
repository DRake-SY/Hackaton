"""
Protocole d'analyse RNAseq :
Ce snakefile prend en entrée une liste de clé d'accession et retourne en sortie une table de comptage ainsi que des fichiers intermédiaires
"""

samples=["SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589"]#  liste des clé d'accession aux données RNAseq


## rule reliant chaque module (rule) entre eux en spécifiant les différents résultats attendue à la fin
rule all:
    input:
        expand(["fastqc/{SAMPLE}_1_fastqc.html","fastqc/{SAMPLE}_2_fastqc.html","chromosome/ref.fa","chromosome/chr_annotation.gtf","star/{SAMPLE}.bam.bai","result.counts"], SAMPLE=samples)
        #on utilise la fonction pour utiliser le wildcard causer par l'utilisation des clés d'accession


## rule permettant de télécharger les données RNAseq utilisées pour l'analyse 
rule download_fastq:
 output:
    "samples/{SAMPLE}_1.fastq","samples/{SAMPLE}_2.fastq"
 singularity:
  "docker://pegi3s/sratoolkit:2.10.0" #docker choisit pour utiliser les commandes de sra-toolkit
 threads: 16
 resources: load = 25  #limte de ressources
 shell: 
  "fasterq-dump {wildcards.SAMPLE} --progress -O samples"



## Rule utilisant la fonction fastqc : elle prend en entrée une paire de fichier fastq et retourne des fichiers Html
rule fastqc:
    input:
        sample1="samples/{SAMPLE}_1.fastq",
        sample2="samples/{SAMPLE}_2.fastq"
    output:
        html1="fastqc/{SAMPLE}_1_fastqc.html",
        html2="fastqc/{SAMPLE}_2_fastqc.html"
    threads: 2
    singularity: #docker créé pour utiliser la fonction fastqc
        "docker://drakesy/hackaton:fastqcv2"
    shell: "fastqc -o fastqc -t {threads} {input.sample1} {input.sample2}"


## Rule qui va télécharger l'ensemble des chromosomes nécessaire à l'indexation star
rule download_chromosome:
   output:
    "chromosome/ref.fa" #fichier regroupant l'ensemble des données de chaque chromosome
   shell:
    """
    for chromosome in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "MT" "X" "Y";
    do
        wget -O "$chromosome".fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome."$chromosome".fa.gz
    done
    gunzip -c *.fa.gz > {output} \
    && rm *.fa.gz
    """

## Rule permettant de télécharge le génome de référence 
rule download_genome_annotation:
    output: "chromosome/chr_annotation.gtf"
    resources: load=25
    shell: "wget -O chromosome/chr_annotation.gtf.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz\
    && gunzip chromosome/chr_annotation.gtf.gz"


## Rule permettant l'indexation du génome de référence 
rule index:
    input:"chromosome/ref.fa" #fichier d'entré des chromosomes
    output:"chromosome/SAindex"
    singularity:"docker://drakesy/hackaton:starv2"
    threads: 16
    resources: load=25
    shell:
     "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir chromosome/ --genomeFastaFiles {input}"


## Rule qui va aligner les reads avec le génome de référence indexé 
rule mapping:
    input: #on utilise la pair de fichier fastq pour chaque clé d'accession
        "chromosome/SAindex",
        sample1="samples/{SAMPLE}_1.fastq", 
        sample2="samples/{SAMPLE}_2.fastq"
    output:"star/{SAMPLE}.bam" #sort des fichiers bam pour chaque clé d'accesion
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

## Rule permettant l'indexation du fichier bam aligné sur le génome de référence
rule samtools:
    input: "star/{SAMPLE}.bam"
    output: "star/{SAMPLE}.bam.bai"
    singularity: "docker://drakesy/hackaton:samtools"
    resources: load=25
    shell: "samtools index {input}"

## Rule permettant la création du table de comptage à partir des fichiers bam
rule counting_reads:
 input:
  samp=expand("star/{smp}.bam", smp=samples),
  gtf="chromosome/chr_annotation.gtf"
 output:
  "result.counts"
 singularity:
  "docker://drakesy/hackaton:subread"
 threads: 4
 resources: load=25
 shell:
  """
   featureCounts -T {threads} -t gene -g gene_id -s 0 -a {input.gtf} -o {output} {input.samp}
  """


