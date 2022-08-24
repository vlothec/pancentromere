# Automated workflow for paired-end ChIP-seq data processing and alignment

This is a Snakemake workflow for automated processing and alignment of paired-end read data derived from chromatin immunoprecipitation followed by high-throughput sequencing (ChIP-seq).

### Requirements

- Installation of [Snakemake](https://snakemake.readthedocs.io/en/stable/) and optionally [conda](https://conda.io/docs/)
- Demultiplexed paired-end reads in gzipped FASTQ format located in the `data/` directory. These should be named according to the following naming convention: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`
- A samtools-indexed reference genome in FASTA format and a chromosome sizes file (e.g., `Col-0.ragtag_scaffolds.fa`, `Col-0.ragtag_scaffolds.fa.fai`, and `Col-0.ragtag_scaffolds.fa.sizes`, the latter two of which generated with `samtools faidx Col-0.ragtag_scaffolds.fa; cut -f1,2 Col-0.ragtag_scaffolds.fa.fai > Col-0.ragtag_scaffolds.fa.sizes`), each located in `data/index/` 
- A reference genome index for bowtie2, located in `data/index/`
- `Snakefile` in this repository. This contains "rules" that each execute a step in the workflow
- `config.yaml` in this repository. This contains customizable parameters including `reference`, which should be the reference genome file name without the `.fa` extension (e.g., `Col-0.ragtag_scaffolds`)
- Optional: `ChIPseq_mapping_environment.yaml` in this repository, used to create the software environment if conda is used
- If conda is not used, the tools listed in ChIPseq_mapping_environment.yaml must be specified in the PATH variable

### Creating the conda environment

```
conda env create --file ChIPseq_mapping_environment.yaml --name ChIPseq_mapping
```

### Usage

In a Unix shell, navigate to the base directory containing `Snakefile`, `config.yaml`, `ChIPseq_mapping_environment.yaml`, and the `data\` subdirectory, which should have a directory tree structure like this:

```
.
├── config.yaml
├── data
│   ├── Col-0_CENH3_PE150_Rep1_ChIP_R1.fastq.gz
│   ├── Col-0_CENH3_PE150_Rep1_ChIP_R2.fastq.gz
│   ├── Col-0_CENH3_PE150_Rep1_input_R1.fastq.gz
│   ├── Col-0_CENH3_PE150_Rep1_input_R2.fastq.gz
│   └── index
│       ├── Col-0.ragtag_scaffolds.1.bt2l
│       ├── Col-0.ragtag_scaffolds.2.bt2l
│       ├── Col-0.ragtag_scaffolds.3.bt2l
│       ├── Col-0.ragtag_scaffolds.4.bt2l
│       ├── Col-0.ragtag_scaffolds.fa
│       ├── Col-0.ragtag_scaffolds.fa.fai
│       ├── Col-0.ragtag_scaffolds.fa.sizes
│       ├── Col-0.ragtag_scaffolds.rev.1.bt2l
│       └── Col-0.ragtag_scaffolds.rev.2.bt2l
├── ChIPseq_mapping_environment.yaml
├── README.md
├── scripts
│   ├── deduplicate_gz_python3.py
│   ├── genomeBin_bedgraphToTSV.R
│   └── keepPaired.py
│   
└── Snakefile
```

Then run the following commands in the base directory (`--cores` should match the `THREADS` parameter in `config.yaml`):

```
conda activate ChIPseq_mapping
snakemake -p --cores 48
conda deactivate
```

### Useful Snakemake parameters

- `--cores` specifies the maximum number of threads
- `-n` performs a dry run
- `-p` prints commands
- `--use-conda`
- `--conda-prefix ~/.myconda`
- `--forcerun calc_coverage` forces rerun of a given rule (e.g., `calc_coverage`)

### Updating the conda environment

```
conda env update --file ChIPseq_mapping_environment.yaml --name ChIPseq_mapping
```
