# USAGE:  # ./script_alignHiFi.sh <WORKDIR> <ACC> <CCS_NAME> <CORES>

# EXAMPLE # ./script_alignHiFi.sh <path_to_working_directory> Rabacal-1 22005 8

# Parse parameteres
WORKDIR=$1
ACC=$2
CCS_NAME=$3
CORES=$4

# Declare other variables
WORKTMP=$PWD
CODE=$PWD

REF=$ACC.Chr_scaffolds.fa
CCS=$CCS_NAME.q20.fastq.gz
ChrM=at_mitochondria.fa
ChrC=at_chloroplast.fa

output=output_dir/$ACC
outputTMP=output_dir/$ACC

# Neccessary software
samtools=<SAMTOOLS_BIN>

# Create necessary directories
mkdir -p $output
mkdir -p $outputTMP

# Start
echo ""
echo "Start script"
date

cd $outputTMP/

echo "Adding ChrM and ChrC to REFerence..."
cat $REF $ChrM $ChrC > $ACC.scaffolds.ChrMChrC.fa

$samtools index $ACC.scaffolds.ChrMChrC.fa

#########################################
################ CCS reads ##############
#########################################

# Activate virtual environment
source activate pbmm2_v1.9.0

echo "Aligning CCS reads..."
pbmm2 align --sort -j $CORES --rg "@RG\tID:$ACC\tSM:$ACC" \
	--log-level DEBUG \
	--preset SUBREAD \
	--min-length 5000 \
	$ACC.scaffolds.ChrMChrC.fa $CCS $outputTMP/$ACC.CCS.bam

echo "Filtering out unmapped reads, secondary and supplementary alignments..."
$samtools view -@ $CORES -b -F 2308 $outputTMP/$ACC.CCS.bam Chr1 Chr2 Chr3 Chr4 Chr5 > $outputTMP/$ACC.CCS.Chr.F2308.bam
$samtools index $outputTMP/$ACC.CCS.Chr.F2308.bam

rm $outputTMP/$ACC.CCS.bam*

conda deactivate

##############################################
################ NucFreq PLOTS ###############
##############################################

cd $outputTMP

grep "$ACC" $CODE/Table_S3.bed > $ACC.Table_S3.bed

# Activate virtual environment with appropriate Python version
source activate python_v3.8

echo "Estimating secondary allele for CENtromeres..."
$CODE/NucFreq-0.1/NucPlot.py --bed $outputTMP/$ACC.Table_S3.bed -y 250 --threads $CORES --minobed 2 --obed $outputTMP/$ACC.CCS.Chr.F2308.CEN.bed $outputTMP/$ACC.CCS.Chr.F2308.bam $outputTMP/$ACC.CCS.Chr.F2308.CEN.png

rsync -av $outputTMP/$ACC.*.png $output/
rsync -av $outputTMP/$ACC.*.bed $output/

conda deactivate

#################################################
############## HetDetection TABLE ###############
#################################################

echo "Generating Heterozygosity tables..."
### CCS ###
Rscript --vanilla $CODE/NucFreq-0.1/HetDetection_modified2.R $outputTMP $ACC.CCS.Chr.F2308.CEN.bed

rsync -av $outputTMP/*.bed* $output/

# End
date
echo "End script"
echo ""

