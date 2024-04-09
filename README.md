# Pipeline for processing PLASMIDSAURUS ONT data

> Álvaro Herrero Reiriz

This pipeline is designed to process ONT data from PLASMIDSAURUS to generate BAM files and sashimi plots. An automated version based on Nextflow is also available [here](https://github.com/a-hr/plasmidsaurus_pipeline).

- [Pipeline for processing PLASMIDSAURUS ONT data](#pipeline-for-processing-plasmidsaurus-ont-data)
  - [Inputs](#inputs)
  - [Requirements](#requirements)
  - [Running the pipeline](#running-the-pipeline)


## Inputs

- FASTQ files from PLASMIDSAURUS: expected to be already demultiplexed and trimmed.
- Reference genome in FASTA format.
- Reference annotation in GTF and BED format (if splice-aware alignment is expected).
- All the inputs required by the ggsashimi module, check the [repository](https://github.com/a-hr/sashimiplots) for more information.


## Requirements

In order to create plots, the following software is required:

- HPC cluster:
	- Singularity
	- Go

- Local machine:
  - Docker

In order to pull the Singularity image:
```bash
module load Singularity Go  # Load the modules if using HPC
singularity pull ggsashimi.sif docker://guigolab/ggsashimi
```
> Note: The image is expected to be in the same directory as the script.

The rest of the requirements are available in the [conda environment file](environment.yml). To create the environment, run:

```bash
conda env create -f environment.yml -p /path/to/env
```

## Running the pipeline

The pipeline is divided into three main steps:

1. **QC**: Run `FASTQC` and `multiqc` on the input FASTQ files.
2. **Alignment**: Align the reads to the reference genome using `minimap2` and generate a BAM file.
3. **Sashimi plots**: Generate sashimi plots using `ggsashimi`. For more information, check the [repository](https://github.com/a-hr/sashimiplots)

Before running the pipeline, make sure to edit the input arguments in the first section of the `run_pipeline_X.sh` script.

**alignment options**

- `fastq_dir`: directory containing the FASTQ files, expected to be gzipped (default `fastqs/`)
- `fa_path`: path to the reference FASTA (default `ref_genome/reference.fa`)
- `gtf_path`: path to the reference GTF (default `ref_genome/reference.gtf`)
- `bed_path`: path to the reference BED (default `ref_genome/reference.bed`)

**sashimi plot options**

- `sashimi_config`: path to the `config.tab` file detailing the evetns to plot (default `inputs/config.tab`)
- `palette_path`: path to the `palette.txt` file (default `inputs/palette.txt`)
- `bam_tsv`: path to the input_bams.tsv file (default `inputs/input_bams.tsv`)
> BAMs will be generated in the `$output_dir/bams` directory, so name your BAMs `$output_dir/bams/FASTQ_NAME.bam`

- `sashimi_min_cov`: the minimum number of reads supporting a junction to be drawn (default `3`)
- `sashimi_agg`: the aggregate function for overlay (default `mean_j`)
- `sashimi_alpha`: the transparency level for the density histogram (default `0.6`)
  
> Additional configuration options can be checked in the [repository](https://github.com/a-hr/sashimiplots)

**output options**
- `output_dir`: path to the output directory (created if it doesn't exist) (default `results/`)

The pipeline can be run using the following command:

Local run:
```bash
bash run_pipeline_local.sh
```

Cluster run:
```bash
sbatch run_pipeline_cluster.sh
```