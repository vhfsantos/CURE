#!/bin/bash

set -e
trap ctrl_c SIGINT
trap 'echo - Ignoring a SIGHUP received..' SIGHUP

error_exit() {
	msg=$1
	>&2 echo -e "\033[1;31mERROR: ${msg}\033[0m"
	usage >&2
	exit 2
}

function ctrl_c() {

	echo "- Canceling all parallel jobs..."
	killall parallel
	echo "- Bye"
	exit 1

}
# Help text for the supported command-line interface.
usage() {
cat <<'HELP'
CURE - estimate trees

Usage:
  estimate-trees.sh --estimated-trees <output-dir> [input options] [options]

Input options (provide at least one):
  --gene-region-out <dir>      CURE GeneRegion output directory
  --uce-region-out <dir>       CURE UCERegion output directory
  --custom-alignments <dir>    Directory of FASTA, PHYLIP, or NEXUS alignments
                               (.fa, .fasta, .phy, .phylip, .nex, .nexus)

Required output option:
  --estimated-trees <dir>      Directory for IQ-TREE results and ASTRAL input

Other options:
  --threads <number>           Maximum concurrent tree estimations (default: 2)
  --only-by-gene               Skip GeneRegion alignments concatenated by genic
                               region; intergenic alignments are still included
  --only-by-genic-region       Skip GeneRegion alignments concatenated by gene;
                               intergenic alignments are still included
  -h, --help                   Show this help and exit

Examples:
  estimate-trees.sh --gene-region-out CURE-GeneRegion-output \
    --uce-region-out CURE-UCERegion-output --estimated-trees estimated-trees

  estimate-trees.sh --gene-region-out CURE-GeneRegion-output \
    --estimated-trees estimated-trees --only-by-gene

  estimate-trees.sh --custom-alignments input-alignments \
    --estimated-trees estimated-trees

The input options can be combined. For CURE inputs, keep the directory and file
names produced by GeneRegion and UCERegion. Requires IQ-TREE and GNU Parallel
in the active CURE environment. The two --only-by-* options require
--gene-region-out and cannot be combined.
HELP
}

# Option strings for arg parser
SHORT=h
LONG=help,gene-region-out:,uce-region-out:,estimated-trees:,threads:,only-by-gene,only-by-genic-region,custom-alignments:

# Set deafult values
tmp=$(realpath "$0")
HOME_DIR=${tmp%/*}
ONLY_BY_GENE="False"
ONLY_BY_REGION="False"
THREADS=2

check_deps() {
	for app in $CONDA_PREFIX/bin/iqtree $CONDA_PREFIX/bin/parallel; do
        	command -v $app >/dev/null 2>&1 || \
			error_exit "Cannot find ${app} in your PATH variable\nDid you activate the CURE environment?"
	done
}


# Read options
if ! OPTS=$(getopt --options "$SHORT" --long "$LONG" --name "$0" -- "$@"); then
	error_exit "Failed to parse options...exiting."
fi
eval set -- "$OPTS"

# Extracting arguments
while true ; do
	case "$1" in
		-h | --help )
		usage
		exit 0
		;;
		--gene-region-out )
		GeneRegion_OUT="$2"
		shift 2
		;;
		--uce-region-out )
		UCERegion_OUT="$2"
		shift 2
		;;
		--estimated-trees )
		OUT="$2"
		shift 2
		;;
		--threads )
		THREADS="$2"
		shift 2
		;;
		--only-by-gene )
		ONLY_BY_GENE="True"
		shift
		;;
		--only-by-genic-region )
		ONLY_BY_REGION="True"
		shift
		;;
		--custom-alignments )
		CUSTOM_ALI="$2"
		shift 2
		;;
		-- )
		shift
		break
		;;
		*)
		error_exit "Please, supply all arguments correctly."
		;;
	esac
done

check_deps

# Checking if any required args is empty
if [ -z "${OUT}" ]; then
	error_exit "Please, supply the output directory name."
fi

# Check if none input was provided
if [ -z "${GeneRegion_OUT}" ] && [ -z "${UCERegion_OUT}" ] && [ -z "${CUSTOM_ALI}" ]; then
 	error_exit "Please, supply the input alignments directory."
fi

if [ "$ONLY_BY_GENE" == "True" ] && [ "$ONLY_BY_REGION" == "True" ]; then
	error_exit "Choose either --only-by-gene or --only-by-genic-region."
fi
if [ -z "${GeneRegion_OUT}" ] && { [ "$ONLY_BY_GENE" == "True" ] || [ "$ONLY_BY_REGION" == "True" ]; }; then
	error_exit "The --only-by-* options require --gene-region-out."
fi

#creating outputs
IQTREE_OUT="$OUT"/IQtree_output
ASTRAL_IN="$OUT"/ASTRAL_input
mkdir -p "$IQTREE_OUT" "$ASTRAL_IN"
#-------- GeneRegion
GR_BY_GENE="$ASTRAL_IN"/GeneRegion_by-gene
GR_BY_GENIC_REG="$ASTRAL_IN"/GeneRegion_by-genic-region
#-------- UCERegion
UR_ASTRAL="$ASTRAL_IN"/UCERegion

# Func to estimate intergenic regions (GeneRegion)
# Run as: input_dir output_dir task
Run_IQtree_NEXUS() {
	NEXUS_DIR="$1"
	OUT_DIR="$2"
	TODO="$3"
	mkdir -p "${OUT_DIR}"
	NEXUS=$(find "$NEXUS_DIR" -name "*.nexus")
	N_NEXUS=$(find "$NEXUS_DIR" -name "*.nexus" | wc -l)
	AUX=0
	if [ "$TODO" == "Run" ]; then
        	echo "- Estimating trees for files in $NEXUS_DIR"
	        echo "- Saving in $OUT_DIR"
	fi
	# run
	for alignment in $NEXUS; do
	        #echo "running $alignment"
		AUX=$(( $AUX+1 ))
		file=$(basename "$alignment")
		$CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs $THREADS \
			$CONDA_PREFIX/bin/iqtree -s "$alignment" --prefix "$OUT_DIR"/${file%.*} \
			-bb 1000 -alrt 1000 --quiet --threads-max 1
		bash ${HOME_DIR}/progress-bar.sh $AUX $N_NEXUS
	done
	$CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
	echo "- Done"
}

# Func to estimate with other regions (GeneRegion) and partitioned UCEs (UCERegion)
# Run as: input_dir output_dir task
Run_IQtree_PHYLIP() {
	PHYLIP_DIR="$1"
	OUT_DIR="$2"
	TODO="$3"
	mkdir -p "${OUT_DIR}"
	PHYLIP=$(find "$PHYLIP_DIR" -name "*.phylip")
	N_PHYLIP=$(find "$PHYLIP_DIR" -name "*.phylip" | wc -l)
	AUX=0
	if [ "$TODO" == "Run" ]; then
	        echo "- Estimating trees for files in $PHYLIP_DIR"
	        echo "- Saving in $OUT_DIR"
	fi
	# run
	for alignment in $PHYLIP; do
		AUX=$(( $AUX+1 ))
		file=$(basename "$alignment")
		$CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs $THREADS \
			$CONDA_PREFIX/bin/iqtree -s "$alignment" -spp "$PHYLIP_DIR"/${file%.*}.charsets \
			--quiet --prefix "$OUT_DIR"/${file%.*} \
			-bb 1000 --threads-max 1
		bash ${HOME_DIR}/progress-bar.sh $AUX $N_PHYLIP
	done
	$CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
	echo "- Done"
}

Run_IQtree_CUSTOM() {
	NEXUS_DIR="$CUSTOM_ALI"
	OUT_DIR="$IQTREE_OUT"/
	TODO="$1"
	mkdir -p "${OUT_DIR}"
	NEXUS=$(find "$NEXUS_DIR" -name "*.fa" -o -name "*.fasta" -o -name "*.phy"\
	 -o -name "*.nex" -o -name "*.nexus" -o -name "*.phylip")
	N_NEXUS=$(echo $NEXUS | wc -w)
	if [ $N_NEXUS -eq 0 ]; then error_exit "Custom alignments must be Phylip, Nexus or Fasta. Make sure the files have correct extensions."; fi
	AUX=0
	if [ "$TODO" == "Run" ]; then
        	echo "- Estimating trees for files in $NEXUS_DIR"
	        echo "- Saving in $OUT_DIR"
	fi
	# run
	for alignment in $NEXUS; do
		AUX=$(( $AUX+1 ))
		file=$(basename "$alignment")
		$CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs $THREADS \
			$CONDA_PREFIX/bin/iqtree -s "$alignment" --quiet \
			--prefix "$OUT_DIR"/${file%.*} -bb 1000 -alrt 1000 \
			--threads-max 1
		bash ${HOME_DIR}/progress-bar.sh $AUX $N_NEXUS
	done
	$CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
	echo "- Done"
}



# Run IQtree for custom data, if provided
if [ -z "$CUSTOM_ALI" ]; then
	echo "- No custom alignment were provided. Skipping..."
else
    echo "- Running IQ-Tree for custom alignment files"
	echo "--------------------------------------------"
	Run_IQtree_CUSTOM Run
	echo "- Checking..."
	Run_IQtree_CUSTOM Check
	# ------------------------------------------------------
	echo "- Preparing ASTRAL input for custom alignments"
	cat "$IQTREE_OUT"/*treefile \
		> "$ASTRAL_IN"/concatenated-trees.tre
fi

# Run IQtree for GeneRegion data, if provided

if [ -z "$GeneRegion_OUT" ]; then
	echo "- No alignments of CURE GeneRegion were provided. Skipping..."
else
    echo "- Running IQ-Tree for CURE GeneRegion files"
	echo "-------------------------------------------"
	# Running for intergenic regions
	mkdir -p "$IQTREE_OUT"/"intergenic-regions"
	if [ -z "$(ls -A "$IQTREE_OUT"/"intergenic-regions")" ]; then
		Run_IQtree_NEXUS "$GeneRegion_OUT"/intergenic-regions \
		                 "$IQTREE_OUT"/intergenic-regions \
						 Run
		echo "- Checking..."
		Run_IQtree_NEXUS "$GeneRegion_OUT"/intergenic-regions \
		                 "$IQTREE_OUT"/intergenic-regions \
						 Check
	else
		echo "- Already estimated trees of intergenic regions. Skipping..."
	fi

	# Running for concatenated by genic region
	if [ "$ONLY_BY_GENE" == "False" ]; then
	        mkdir -p "$IQTREE_OUT"/"concatenated-by-genic-region"
		if [ -z "$(ls -A "$IQTREE_OUT"/"concatenated-by-genic-region")" ]; then
			Run_IQtree_NEXUS "$GeneRegion_OUT"/concatenated-by-genic-region \
			                 "$IQTREE_OUT"/concatenated-by-genic-region \
							 Run
			echo "- Checking..."
			Run_IQtree_NEXUS "$GeneRegion_OUT"/concatenated-by-genic-region \
			                 "$IQTREE_OUT"/concatenated-by-genic-region \
							 Check
			# ------------------------------------------------------
			echo "- Preparing ASTRAL input (by genic region)"
			mkdir -p $GR_BY_GENIC_REG
			# only exons
			cat "$IQTREE_OUT"/concatenated-by-genic-region/*__exon*treefile \
				> "$GR_BY_GENIC_REG"/only-exons.tre
			# only introns
			cat "$IQTREE_OUT"/concatenated-by-genic-region/*__introns.treefile \
				> "$GR_BY_GENIC_REG"/only-introns.tre
			# only genic regions (exons + introns)
			cat "$GR_BY_GENIC_REG"/only-introns.tre "$GR_BY_GENIC_REG"/only-exons.tre \
				> "$GR_BY_GENIC_REG"/only-exons-and-introns.tre
			# all regions (exons + introns + intergenic regions)
			cat "$IQTREE_OUT"/intergenic-regions/*treefile \
				"$GR_BY_GENIC_REG"/only-introns.tre "$GR_BY_GENIC_REG"/only-exons.tre \
				> "$GR_BY_GENIC_REG"/all-regions.tre

		else
			echo "- Already estimated trees by genic region. Skipping..."
		fi
	fi

	# Create output concatenated by gene
	if [ "$ONLY_BY_REGION" == "False" ]; then
	        mkdir -p "$IQTREE_OUT"/"concatenated-by-gene"
		if [ -z "$(ls -A "$IQTREE_OUT"/"concatenated-by-gene")" ]; then
		    # $IQTREE_OUT"/"
			Run_IQtree_PHYLIP "$GeneRegion_OUT"/concatenated-by-gene \
			                  "$IQTREE_OUT"/concatenated-by-gene \
							  Run
			echo "- Checking..."
			Run_IQtree_PHYLIP "$GeneRegion_OUT"/concatenated-by-gene \
			                  "$IQTREE_OUT"/concatenated-by-gene \
							  Check
			# ------------------------------------------------------
			echo "- Preparing ASTRAL input (by gene)"
			mkdir -p $GR_BY_GENE
			# only genic regions
			cat "$IQTREE_OUT"/concatenated-by-gene/*treefile \
				> "$GR_BY_GENE"/only-exons-and-introns.tre
			# all regions (add intergenic)
			cat "$IQTREE_OUT"/intergenic-regions/*treefile \
				"$GR_BY_GENE"/only-exons-and-introns.tre \
				> "$GR_BY_GENE"/all-regions.tre
		else
			echo "- Already estimated trees by gene. Skipping..."
		fi
	fi
fi

# Run IQtree for UCERegion data, if provided

if [ -z "$UCERegion_OUT" ]; then
	echo "- No alignments of CURE UCERegion were provided. Skipping..."
else
    echo "- Running IQ-Tree for CURE UCERegion files"
	echo "------------------------------------------"
	
	Run_IQtree_PHYLIP "$UCERegion_OUT"/partitioned-uces \
		             "$IQTREE_OUT"/UCERegion \
					 Run
	echo "- Checking..."
	Run_IQtree_PHYLIP "$UCERegion_OUT"/partitioned-uces \
		             "$IQTREE_OUT"/UCERegion \
					 Check
	# ------------------------------------------------------
	mkdir -p "$UR_ASTRAL"
	echo "- Preparing ASTRAL input for UCERegion files"
	cat "$IQTREE_OUT"/UCERegion/*treefile \
		> "$UR_ASTRAL"/partitioned-uces.tre

fi
echo "- All done!

If you make use of CURE in your research, please cite:

Freitas et al. (2023). UCE phylogenomics, biogeography, and classification of long-horned
bees (Hymenoptera: Apidae: Eucerini), with insights on using specimens with extremely 
degraded DNA. Insect Systematics and Diversity 7(4):3.
https://doi.org/10.1093/isd/ixad012"
