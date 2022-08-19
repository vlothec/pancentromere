# USAGE:  # ./script_multiple_alignment.sh <WORKDIR> <SC> <gene_list> <CORES>

# EXAMPLE # ./script_multiple_alignment.sh <path_to_working_directory> 0.98 gene_list01.txt 10 

# Parse parameteres
WORKDIR=$1
SC=$2
gene_list=$3
CORES=$4

# Declare other variables
WORKTMP=$PWD

# Software required
samtools=<path_to_SAMTOOLS_bin>

input=$PWD
inputTMP=$PWD

# Start
echo ""
echo "Start script"
date

# Loop through the list of GENES line by line
while IFS= read -r line; do

	GENE=$(echo "${line}" | cut -f1)
	NAME=$(echo "${line}" | cut -f2)

	gID="$GENE"_"$NAME"

	output=$PWD/$ACC/$gID
	outputTMP=$PWD/$ACC/$gID

	# Create necessary directories
	mkdir -p $output
	mkdir -p $outputTMP

	cd $outputTMP/

	echo "Detecting extra copies for $gID..." 
	grep "$GENE" $inputTMP/*/*.genes.gene.gff | grep -v extra_copy_number=0 > $outputTMP/$gID.extra_copy_number.txt
	cat $outputTMP/$gID.extra_copy_number.txt

	echo "Extracting $gID for each accession..."

	COUNTER=0
	for i in $(ls $inputTMP/*/*.CDS.fa); do
		COUNTER=$[$COUNTER +1]

		ACC=$(echo $i | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
		
		$samtools faidx $inputTMP/$ACC/$ACC.CDS.fa $GENE.1 > $outputTMP/$GENE.$ACC.nt.fa

		sed -i "s/$GENE.1/$ACC/" $outputTMP/$GENE.$ACC.nt.fa
	done

	echo "Concatenating $gID..."
	cat $outputTMP/$GENE.*.nt.fa > $outputTMP/$GENE."$COUNTER"accs.nt.fa

	echo "Running multiple alignment for $gID..."

	# Activate virtual environment
	source activate mafft

	mafft --thread 5 --threadtb 10 --threadit 10 --reorder --maxiterate 1000 --retree 1 --genafpair $outputTMP/$GENE."$COUNTER"accs.nt.fa > $outputTMP/$GENE."$COUNTER"accs.nt.mafft.fa

	conda deactivate

        echo "Transform Fastas into Nexus"
	# Activate virtual environment
	source activate pgdspider_v2.1.1.5

	PGDSpider2-cli -inputfile $outputTMP/$GENE."$COUNTER"accs.nt.mafft.fa -inputformat FASTA -outputfile $outputTMP/$GENE."$COUNTER"accs.nt.mafft.nex -outputformat NEXUS -spid $WORKTMP/Multiple_alignment.nt.spid
	
	conda deactivate

	# Transfer relevant files to server
	rsync -av $outputTMP/$GENE."$COUNTER"accs.nt.fa $output/
	rsync -av $outputTMP/$GENE."$COUNTER"accs.nt.mafft.fa $output/

done < $WORKDIR/$gene_list

# End
date
echo "End script"
echo ""

