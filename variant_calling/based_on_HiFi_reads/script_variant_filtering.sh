#!/bin/bash

# Define variables
REF=Col-0.Chr1to5.ChrMC.fa
REPEATS=Col-0.Repeats_merged.gff

# Load environment and software
source activate bcftools_v1.15.1
GATK4=<path_to_gatk-4.1.3.0>/gatk

# Add all possible INFO tags
bcftools +fill-tags 66_HiFi_genomes.bcf -Ob -o 66_HiFi_genomes.INFO.bcf -- -t all,AN,F_MISSING,TYPE,'INFO/COV=sum(FORMAT/DP)','INFO/AD=sum(FORMAT/AD)'

# Convert BCF to VCF
bcftools view 66_HiFi_genomes.INFO.bcf | bgzip -@ 4 -c > 66_HiFi_genomes.INFO.vcf.gz

#index VCF file
tabix -p vcf 66_HiFi_genomes.INFO.vcf.gz

# Count number of variants left at this point
zcat 66_HiFi_genomes.INFO.vcf.gz | grep -c -v "^#"
# 8213281

# Represent repeats as intervals
awk '{OFS="\t";print $1,$4,$5}' $REPEATS | grep -v "#" | grep "^Chr" > repeats.bed

# Filter out variants overlapping intervals (which also excludes organellar chromosomes)
$GATK4 --java-options '-Xmx16G -Xms16G' VariantFiltration \
  -R $REF \
  -V 66_HiFi_genomes.INFO.vcf.gz \
  -O 66_HiFi_genomes.INFO.minusREPEATS.vcf.gz \
  --exclude-intervals repeats.bed

# Count number of variants left at this point
zcat 66_HiFi_genomes.INFO.minusREPEATS.vcf.gz | grep -c -v "^#"
# 7359371

# Create table from INFO field for further inspection in R
bcftools query -f '%CHROM %POS %REF %ALT %QUAL %TYPE %NS %AN %F_MISSING %COV %AD %AC %AF %AC_Hom %AC_Het %ExcHet %HWE\n' 66_HiFi_genomes.INFO.minusREPEATS.vcf.gz > 66_HiFi_genomes.INFO.minusREPEATS.txt
wc -l 66_HiFi_genomes.INFO.minusREPEATS.txt
# 7359371

###
### At this point, I took it to R to determine the filtering parameters
###


# Keep only biallelic SNPs and INDELs
bcftools view -O z -o 66_HiFi_genomes.INFO.minusREPEATS.biallelicSNPsINDELs.vcf.gz -m2 -M2 66_HiFi_genomes.INFO.minusREPEATS.vcf.gz   
zcat 66_HiFi_genomes.INFO.minusREPEATS.biallelicSNPsINDELs.vcf.gz | grep -c "^Chr"
# 6790003

# Exclude based on ...
bcftools view -O z -o 66_HiFi_genomes.INFO.minusREPEATS.biallelicSNPsINDELs.filter95.vcf.gz -e 'NS <= 62 || F_MISSING > 0.05 || AC = 1 || COV < 2604 || COV > 6685 || AC_Het >= 2' 66_HiFi_genomes.INFO.minusREPEATS.biallelicSNPsINDELs.vcf.gz
zcat 66_HiFi_genomes.INFO.minusREPEATS.biallelicSNPsINDELs.filter95.vcf.gz | grep -c "^Chr"
# 5194221

### OR, if we want to be strciter ###

# Exclude based on ...
bcftools view -O z -o 66_HiFi_genomes.INFO.minusREPEATS.biallelicSNPsINDELs.filter100.vcf.gz -e 'NS < 66 || F_MISSING > 0 || AC = 1 || COV < 2604 || COV > 6685 || AC_Het >= 2' 66_HiFi_genomes.INFO.minusREPEATS.biallelicSNPsINDELs.vcf.gz
zcat 66_HiFi_genomes.INFO.minusREPEATS.biallelicSNPsINDELs.filter100.vcf.gz | grep -c "^Chr"
# 4436685

#Generate useful tables for analysis in R
bcftools query -f '%CHROM %POS[\t%GT]\n' 66_HiFi_genomes.INFO.minusREPEATS.biallelicSNPsINDELs.filter95.vcf.gz > 66_HiFi_genomes.INFO.minusREPEATS.biallelicSNPsINDELs.filter95.GENOTYPES.txt

bcftools query -f '%CHROM %POS[\t%GT]\n' 66_HiFi_genomes.INFO.minusREPEATS.biallelicSNPsINDELs.filter100.vcf.gz > 66_HiFi_genomes.INFO.minusREPEATS.biallelicSNPsINDELs.filter100.GENOTYPES.txt


