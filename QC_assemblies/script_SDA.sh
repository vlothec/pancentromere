# USAGE:  # ./script_SDA.sh <WORKDIR> <ACC> <CORES>

# EXAMPLE # ./script_SDA.sh <path_to_working_directory> Rabacal-1 24

# Parse parameteres
WORKDIR=$1
ACC=$2
CORES=$3

# Declare other variables
WORKTMP=$PWD
CODE=$PWD
SCRIPTS=<path_to_SDA_environment>/sandbox/SDA/scripts

REF=$ACC.Chr_scaffolds.fa

output=output_dir/$ACC
outputTMP=output_dir/$ACC

# Neccessary software
samtools=<SAMTOOLS_BIN>
BEDTOOLS=<BEDTOOLS_BIN>

# Create necessary directories
mkdir -p $output
mkdir -p $outputTMP

# Start
echo ""
echo "Start script"
date

cd $outputTMP/

echo "Soft links with generic names to input files..." 
ln -s $ACC.CCS.Chr.F2308.bam input.bam
ln -s $ACC.CCS.Chr.F2308.bam.bai input.bam.bai
ln -s $REF reference.fa
$samtools faidx reference.fa

#############################################
################## SDA cluster ##############
#############################################

echo "Removing previous outputs and snakemake..."
rm -rf sda_out
rm -rf .snakemake
rm -rf RM_*

### Activate virtual environment
source activate sda_4ec1b3a

echo "Executing single command..."
singularity exec ${SINGULARITY_SANDBOX} SDA denovo \
	--threads 24 \
	--platform subread \
	--species library.rDNA_telomeres.fa \
	--input input.bam \
	--ref reference.fa \
	--pre sda

echo "Check that we are still in the SDA environment"
conda deactivate

#### BEDtools ####

cd $outputTMP/sda_out/coverage

grep "$ACC" $CODE/Table_S3.bed > $ACC.Table_S3.bed

$BEDTOOLS intersect -a sda.collapses.with.cm.bed -b $ACC.Table_S3.bed -wa > sda.collapses.with.cm.CEN.bed

# Virtual environmet with the right python
source activate python_v3.8

mkdir -p $output/sda_out/coverage
cd $outputTMP/sda_out/coverage

# Extract expected coverage
COV=$(awk 'FNR == 2 {print $1}' sda.coverage.stats)
echo $COV

$SCRIPTS/count_collapse.py -c $COV sda.collapses.with.cm.CEN.bed > sda.collapses.with.cm.CEN.txt

rsync -av * $output/sda_out/coverage/

conda deactivate

# End
date
echo "End script"
echo ""

