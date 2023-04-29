#!/bin/bash

source activate glnexus_v1.4.1

WORKDIR=<path_to_working_directory>

mkdir $WORKDIR/3_joint_genotyping

cd $WORKDIR

mkdir joint_genotyping
cd joint_genotyping

# Run Glnexus
glnexus_cli --config DeepVariant ../variant_calling/*.g.vcf > 66_HiFi_genomes.bcf

# Remove directory with external sorting of the gVCF data
rm -rf GLnexus.DB

bcftools view 66_HiFi_genomes.bcf | bgzip -@ 4 -c > 66_HiFi_genomes.vcf.gz

conda deactivate

