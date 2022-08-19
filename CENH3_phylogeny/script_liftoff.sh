# USAGE:  # ./script_liftoff.sh <WORKDIR> <TARGET> <ACC> <SC> <CORES>

# EXAMPLE # ./script_liftoff.sh <path_to_working_directory> 22005.Chr_scaffolds.fa Rabacal-1 0.98 8

# Parse parameteres
WORKDIR=$1
TARGET=$2
ACC=$3
SC=$4
CORES=$5

# Declare other variables
WORKTMP=$PWD

TAIR_GENES=Araport11_GFF3_genes_transposons.201606.gff

output=$PWD/$ACC
outputTMP=$PWD/$ACC

# Create necessary directories
mkdir -p $output
mkdir -p $outputTMP

# Activate virtual environment
source activate liftoff_v1.6.2

# Start
echo ""
echo "Start script"
date

cd $outputTMP/

echo "Soft link to gff"
ln -s $TAIR_GENES ./TAIR_GENES.gff

echo "Running liftoff annotation ..."

liftoff -g TAIR_GENES.gff \
	-p $CORES \
	-copies -sc $SC \
	-polish -cds \
	-o $outputTMP/$ACC.genes.gff \
	$TARGET $TAIR

mv $outputTMP/unmapped_features.txt $outputTMP/$ACC.unmapped_features.txt

grep -P "\tgene" $outputTMP/$ACC.genes.gff_polished > $outputTMP/$ACC.genes.gene.gff 

conda deactivate

echo "Extracting CDS per gene..."

# Activate virtual environment
source activate gffread_v0.12.7

gffread -x $outputTMP/$ACC.CDS.fa \
	-g $TARGET \
	$outputTMP/$ACC.genes.gff_polished

conda deactivate

# Remove unnecessary files
rm $outputTMP/TAIR_GENES.gff_db
rm $outputTMP/intermediate_files/*

# Transfer to project directory 
rsync -av $outputTMP/$ACC.unmapped_features.txt $output/
rsync -av $outputTMP/$ACC.genes.gff_polished $output/
rsync -av $outputTMP/$ACC.genes.gene.gff $output/
rsync -av $outputTMP/$ACC.*.fa $output/

# End
date
echo "End script"
echo ""

