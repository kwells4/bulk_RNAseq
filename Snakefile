from snakemake.utils import min_version
import pandas as pd

##### set minimum snakemake version #####
min_version("5.1.4")

##### load config file #####
configfile: "config.yaml"
samples = pd.read_table("samples.tsv").set_index("sample", drop = False)
SAMPLE_LIST = samples.index.values
data_dir = config['data_dir']
is_paired = "fastq2" in samples.columns

def get_input(wildcards):
    fastq1 = samples.loc[wildcards.sample, "fastq1"]
    fastq1 = data_dir + fastq1
    if is_paired:
        fastq2 = samples.loc[wildcards.sample, "fastq2"]
        fastq2 = data_dir + fastq2
        return(fastq1, fastq2)
    else:
        return(fastq1)

rule all:
    input:
         expand("aligned/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLE_LIST),
         expand("fastqc/{sample}_fastqc.txt", sample=SAMPLE_LIST)

rule count_all:
    input:
        expand("featureCount/{sample}_countsOutput", sample=SAMPLE_LIST)

rule fastqc:
    input:
        input_list=get_input
    output:
        "fastqc/{sample}_fastqc.txt"
    params:
        output_dir="fastqc/",
        sample_name="{sample}"
    conda:
        "envs/alignment_counting.yaml"
    shell:
        """
        fastqc {input} --outdir {params.output_dir}
        echo {params.sample_name} > {output}
        """    

rule align:
    input:
        input_list=get_input
    output:
        "aligned/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        genome=config['reference'],
        gtf=config['gtf'],
        output_prefix="aligned/{sample}_"
    threads: 10
    conda:
        "envs/alignment_counting.yaml"
    shell:
        "STAR --runThreadN {threads} --genomeDir {params.genome} --sjdbGTFfile "
        "{params.gtf} --readFilesIn {input} --readFilesCommand zcat --outSAMtype "
        "BAM SortedByCoordinate --outFileNamePrefix {params.output_prefix}"

rule count:
    input:
        "aligned/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "featureCount/{sample}_countsOutput"
    params:
        gtf=config['gtf'],
        output_file="countsOutput"
    conda:
        "envs/alignment_counting.yaml"
    shell:
        "featureCounts -a {params.gtf} -o {output} -s 0 -R CORE --primary {input}"

rule count_table:
    input:
        expand("featureCount/{sample}_countsOutput", sample=SAMPLE_LIST)
    output:
        "analysis_outs/{project}_countTable.txt"
    script:
        "scripts/countTable.py"

