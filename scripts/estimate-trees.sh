#!/bin/bash

set -e

# usage
usage() {
echo -e "
-----------------------------------------------------------
CURE: an automated and parallel pipeline for UCE curation
-----------------------------------------------------------
--------------------- estimate trees ----------------------
-----------------------------------------------------------
by Vinícius H F Santos & Felipe V Freitas

\e[4mUsage\e[0m: 
 estimate-trees.sh --input-dir <path/to/input-dir> \\
                   --output-dir <path/to/output-dir>

\e[4mRequired arguments\e[0m:
 --input-dir            Path to alignments produced by CURE
 --output-dir           Output directory name

\e[4mOptional arguments\e[0m:
 --threads              Number of threads for the analysis (Default: 10)
 --only-by-gene         Only estimate trees for 'concatenated-by-gene/' files
 --only-by-region       Only estimate trees for 'concatenated-by-region/' files
 

 ╔══════════════════════════════════════════════════════╗
 ║                          NOTE                        ║
 ╠══════════════════════════════════════════════════════╣ 
 ║  The parameter 'input-dir' must direct to the output ║
 ║  of CURE. This script will enter this directory and  ║
 ║  look for the outputs produced by CURE.              ║
 ║  If you renamed them, it will raise an error.        ║
 ╚══════════════════════════════════════════════════════╝
 "
exit 2
}

# Option strings for arg parser
SHORT=h
LONG=help,input-dir:,output-dir:,threads:,only-by-gene,only-by-region

# Set deafult values
tmp=$(realpath "$0")
HOME_DIR=${tmp%/*}
ONLY_BY_GENE="False"
ONLY_BY_REGION="False"
THREADS=10

check_deps() {
	for app in iqtree parallel; do
        	command -v $app >/dev/null 2>&1 || \
			error_exit "Cannot find ${app} in your PATH variable\nDid you activate the cure environment?"
	done
}

# Read options
OPTS=$(getopt --options $SHORT --long $LONG --name "$0" -- "$@")
if [ $? != 0 ]; then error_exit "Failed to parse options...exiting."; fi
eval set -- "$OPTS"

# Extracting arguments
while true ; do
	case "$1" in
		-h | --help )
		usage
		;;
		--input-dir )
		INPUT_DIR="$2"
		shift 2
		;;
		--output-dir )
		OUTPUT_DIR="$2"
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

Run_IQtree_NEXUS() {
	NEXUS_DIR="$INPUT_DIR"/"$1"
	OUTPUT_DIR="$OUTPUT_DIR"/"$1"
	mkdir -p "${OUTPUT_DIR}"
	NEXUS=$(find "$NEXUS_DIR" -name *.nexus)
	N_NEXUS=$(find "$NEXUS_DIR" -name *.nexus | wc -l)
	AUX=0
	echo "$NEXUS_DIR"
	echo "$OUTPUT_DIR"

#	# run
#for alignment in $NEXUS; do
#		file=$(basename "$alignment")
#		iqtree -s "$alignment" --quiet --prefix
}

Run_IQtree_NEXUS intergenic-regions