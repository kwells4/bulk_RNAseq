This snakemake pipeline contains all the necessary scripts to recreate the analysis for (eventual paper).

Writen by Kristen Wells July 16, 2019. Last modified July 16, 2019

To use:

1. Download and install miniconda3: For Linux

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh bash Miniconda3-latest-Linux-x86_64.sh
```

2. Install Snakemake:

```bash
conda install snakemake -c bioconda -c conda-forge
```

3. Update the config file (config.yaml)

* reference: The path to the directory created using STAR
* GTF: The path to the gtf file. Should match the genome fa
* data_dir: Location of the raw data

4. Update the samples.tsv file to include the samples and file names of the samples. This document should have a column named "sample" and a column named "fastq1". If you have paired end reads, also include a column named "fastq2". This file must be
tab delimited.

5. If you are running on a cluster, update the cluster.json for with your specific specs

The command to align all files is: snakemake

Nothing else if required to align all files.

The command to align one file is: snakemake {sample_name}_Aligned.sortedByCoord.out.bam

I recommend submitting to the cluster. Example submit scripts are "alignment.sh" (to align all samples in the sample sheet) or "align_one.sh" to align just one. These are both scripts that will submit jobs to the cluster. The final line is the submit command. This command tells snakemake how to submit jobs to the cluster. It also includes the command to perform the alignment.

