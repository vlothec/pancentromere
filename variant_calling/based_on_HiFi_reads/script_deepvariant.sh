# USAGE:  # ./script_deepvariant.sh <WORKDIR> <ACC> <CORES>
# EXAMPLE # ./script_deepvariant.sh <path_to_working_directory> Rabacal-1 8

# Parse parameteres
WORKDIR=$1
ACC=$2
CORES=$3

# Declare other variables
output=$WORKDIR/variant_calling

# Create necessary directories
mkdir -p $output

# Start
echo ""
echo "Start script"
date

cd $output/

#######################################################
############## DeepVariants with udocker ##############
#######################################################

# Activate virtual environment
source activate udocker

# Add this line when running on cluster!
export PROOT_NO_SECCOMP=1
export OMP_NUM_THREADS=1

INPUT_DIR=$WORKDIR
OUTPUT_DIR=$WORKDIR

CONTAINER_ID=dv_$ACC
ref=Col-0

echo "Creating container..."

udocker create --name=$CONTAINER_ID google/deepvariant:1.3.0

echo "Running DeepVariant..."

udocker run \
	-v "${INPUT_DIR}":"/input" \
	-v "${OUTPUT_DIR}":"/output" \
	$CONTAINER_ID \
	/opt/deepvariant/bin/run_deepvariant \
	--model_type=PACBIO \
	--ref=/input/ref/$ref.Chr1to5.ChrMC.fa \
	--reads=/input/remapping_CCS/$ACC/$ACC.$ref.CCS.bam \
	--num_shards="$CORES" \
	--output_vcf=/output/variant_calling/$ACC.$ref.DV.vcf \
	--output_gvcf=/output/variant_calling/$ACC.$ref.DV.g.vcf 

date

sleep 10

echo "Uninstalling deepvariant specific container..."

cd /tmp/.udocker/containers/
udocker rm $CONTAINER_ID

conda deactivate


# End
date
echo "End script"
echo ""

