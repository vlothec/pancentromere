# Coverage profiles around genomic features

This subdirectoy contains a Snakemake workflow for creating a matrix of windowed coverage values within genomic features and in flanking regions.


### Requirements

- Installation of [Snakemake](https://snakemake.readthedocs.io/en/stable/) and optionally [conda](https://conda.io/docs/)
- `Snakefile` in this repository. This contains "rules" that each execute a step in the workflow
- `config.yaml` in this repository. This contains customizable parameters including `reference`, which should be the reference genome file name without the `.fa` extension (e.g., `Col-0.ragtag_scaffolds`)
- Optional: `ChIPseq_mapping_environment.yaml`, used to create the software environment if conda is used
- If conda is not used, [deepTools](https://deeptools.readthedocs.io/en/develop/) must be installed and specified in the PATH variable
- Genomic feature coordinates and, separately, random locus coordinates in BED6 format: column 1 = chromosome ID; column 2 = 0-based start coordinates; column 3 = 1-based end coordinates; column 4 = sequential or otherwise unique numbers (this speeds up computation; see comment from dpryan79 on 13/09/2018 under GitHub issue [computeMatrix has problem with multi processors #760](https://github.com/deeptools/deepTools/issues/760)]); column 5 = fill with NA; column 6 = fill with \*
- A bigWig coverage file (generated using deepTools bamCoverage as part of the read alignment pipeline) to be used for calculating coverage profiles around genomic features and random loci (e.g., `Col-0_CENH3_PE150_Rep1_ChIP_MappedOn_Col-0.ragtag_scaffolds_lowXM_both_sort_norm.bw`)

### Creating the conda environment

```
conda env create --file environment.yaml --name ChIPseq_mapping
```

### Usage

In a Unix shell, navigate to the base directory containing `Snakefile` and  `config.yaml`.
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
conda env update --file environment.yaml --name ChIPseq_mapping
```
