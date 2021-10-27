#!/bin/bash

set -e
trap ctrl_c SIGINT
trap 'echo - Ignoring a SIGHUP received..' SIGHUP

error_exit() {
	msg=$1
	>&2 echo -e "\033[1;31mERROR: ${msg}\033[0m"
	usage
}

function ctrl_c() {

	echo "- Canceling all parallel jobs..."
	killall parallel
	echo "- Bye"
	exit 1

}
# usage
usage() {
echo -e "
-----------------------------------------------------------
CURE: an automated and parallel pipeline for UCE curation
-----------------------------------------------------------
--------------------- estimate trees ----------------------
-----------------------------------------------------------
by Vinícius Franceshini-Santos & Felipe Freitas

\e[4mUsage\e[0m:
 estimate-trees.sh --cure-out <path/to/cure-out> \\
                   --iqtree-out <path/to/iqtree-out>

\e[4mRequired arguments\e[0m:
 --cure-out             Path to alignments produced by CURE
 --iqtree-out           Output directory name

\e[4mOptional arguments\e[0m:
 --threads              Number of threads for the analysis (Default: 2)
 --only-by-gene         Only estimate trees for 'concatenated-by-gene/' files
 --only-by-region       Only estimate trees for 'concatenated-by-region/' files
 --custom-alignments    Path to a custom alignment dir. Estimate trees using
                        alignments in this dir, ignoring other inputs.


\e[4mExamples\e[0m:

To estimate trees from alignments produced by CURE, run:

    estimate-trees.sh \\
            --cure-out CURE-output/ \\
	    --iqtree-out estimated-trees

To estimate trees only from alignments concatenated by gene, run:

	estimate-trees.sh \\
	    --cure-out CURE-output/ \\
	    --iqtree-out estimated-trees \\
	    --only-by-gene

To estimate trees only from alignments concatenated by region, run:

	estimate-trees.sh \\
	    --cure-out CURE-output/ \\
	    --iqtree-out estimated-trees \\
	    --only-by-region

Also, you can use this script to estimate gene trees in parallel from alignments
produced by any other tool rather than CURE. For so, just use --custom-alignments
parameter:

        estimate-trees.sh \\
	    --custom-alignments input-alignments/ \\
	    --iqtree-out estimated-trees


 ╔══════════════════════════════════════════════════════╗
 ║                          NOTE                        ║
 ╠══════════════════════════════════════════════════════╣
 ║  The parameter 'cure-out' must direct to the output  ║
 ║  of CURE. This script will enter this directory and  ║
 ║  look for the outputs produced by CURE.              ║
 ║  If you renamed them, it will raise an error.        ║
 ╚══════════════════════════════════════════════════════╝
 "
exit 2
}

# Option strings for arg parser
SHORT=h
LONG=help,cure-out:,iqtree-out:,threads:,only-by-gene,only-by-region,custom-alignments:

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
OPTS=$(getopt --options $SHORT --long $LONG --name "$0" -- "$@")
if [ $? != 0 ]; then error_exit "Failed to parse options...exiting."; fi
eval set -- "$OPTS"

if [ $? != 0 ]; then error_exit "Failed to parse options...exiting."; fi

# Extracting arguments
while true ; do
	case "$1" in
		-h | --help )
		usage
		;;
		--cure-out )
		CURE_OUT="$2"
		shift 2
		;;
		--iqtree-out )
		IQTREE_OUT="$2"
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
		--only-by-region )
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
if [ -z "${IQTREE_OUT}" ]; then
	error_exit "Please, supply the output directory name."
fi

if [ -z "${CURE_OUT}" ] && [ -z "${CUSTOM_ALI}" ]; then
 	error_exit "Please, supply the input aligments directory."
fi

Run_IQtree_NEXUS() {
	NEXUS_DIR="$CURE_OUT"/"$1"
	OUT_DIR="$IQTREE_OUT"/"$1"
	TODO="$2"
	mkdir -p "${OUT_DIR}"
	NEXUS=$(find "$NEXUS_DIR" -name *.nexus)
	N_NEXUS=$(find "$NEXUS_DIR" -name *.nexus | wc -l)
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
			$CONDA_PREFIX/bin/iqtree -s "$alignment" --quiet --prefix "$OUT_DIR"/${file%.*} \
			-bb 1000 -alrt 1000 --threads-max 1
		${HOME_DIR}/progress-bar.sh $AUX $N_NEXUS
	done
	$CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
	echo "- Done"
}

Run_IQtree_PHYLIP() {
	PHYLIP_DIR="$CURE_OUT"/"$1"
	OUT_DIR="$IQTREE_OUT"/"$1"
	TODO="$2"
	mkdir -p "${OUT_DIR}"
	PHYLIP=$(find "$PHYLIP_DIR" -name *.phylip)
	N_PHYLIP=$(find "$PHYLIP_DIR" -name *.phylip | wc -l)
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
		${HOME_DIR}/progress-bar.sh $AUX $N_NEXUS
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
		${HOME_DIR}/progress-bar.sh $AUX $N_NEXUS
	done
	$CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
	echo "- Done"
}

if [ -z "$CUSTOM_ALI" ]; then
	# Running for intergenic regions
	if [ -z "$(ls -A "$IQTREE_OUT"/"intergenic-regions")" ]; then
		Run_IQtree_NEXUS intergenic-regions Run
		echo "- Checking..."
		Run_IQtree_NEXUS intergenic-regions Check
	else
		echo "- Already estimated trees of intergenic regions. Skipping..."
	fi

	# Running for concatenated by region
	if [ "$ONLY_BY_GENE" == "False" ]; then
		if [ -z "$(ls -A "$IQTREE_OUT"/"concatenated-by-region")" ]; then
			Run_IQtree_NEXUS concatenated-by-region Run
			echo "- Checking..."
			Run_IQtree_NEXUS concatenated-by-region Check

		else
			echo "- Already estimated trees by region. Skipping..."
		fi
	fi

	# Create output concatenated by gene
	if [ "$ONLY_BY_REGION" == "False" ]; then
		if [ -z "$(ls -A "$IQTREE_OUT"/"concatenated-by-gene")" ]; then
			Run_IQtree_PHYLIP concatenated-by-gene Run
			echo "- Checking..."
			Run_IQtree_PHYLIP concatenated-by-gene Check
		else
			echo "- Already estimated trees by gene. Skipping..."
		fi
	fi
else
	Run_IQtree_CUSTOM Run
	echo "- Checking..."
	Run_IQtree_CUSTOM Check
fi

echo "- All done. Bye"
