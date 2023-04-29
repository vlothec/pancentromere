# USAGE:  # ./script_alignHiFi.sh <WORKDIR> <ACC> <CCS_NAME> <CORES>

# EXAMPLE # ./script_alignHiFi.sh <path_to_working_directory> Rabacal-1 22005 8

# Parse parameteres
WORKDIR=$1
ACC=$2
CCS_NAME=$3
CORES=$4

# Declare other variables

REF=Col-0.Chr1to5.ChrMC.mmi
CCS=lima/$CCS_NAME.q20.fastq.gz

output=$WORKDIR/remapping_CCS/$ACC

# Create necessary directories
mkdir -p $output

# Start
echo ""
echo "Start script"
date

cd $outputTMP/

#######################################
############## CCS reads ##############
#######################################

# Activate virtual environment
source activate pbmm2_v1.9.0

echo "Indexing of $REF already took place." 
##RAN MANUALLY ONCE## pbmm2 index --preset CCS Col-0.Chr1to5.ChrMC.fa Col-0.Chr1to5.ChrMC.mmi

echo "Aligning CCS reads..."
pbmm2 align --sort -j $CORES --rg "@RG\tID:$ACC\tSM:$ACC" \
	--preset CCS \
	$REF $CCS $outputTMP/$ACC.Col-0.CCS.bam

conda deactivate

# End
date
echo "End script"
echo ""

