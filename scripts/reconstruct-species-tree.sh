#!/bin/bash

set -e

error_exit() {
	msg=$1
	>&2 echo -e "\033[1;31mERROR: ${msg}\033[0m"
	usage
}

# usage
usage() {
echo -e "
-----------------------------------------------------------
 CURE: an automated and parallel pipeline for UCE curation
-----------------------------------------------------------
---------------- reconstruct species tree -----------------
-----------------------------------------------------------
by Vinícius H F Santos & Felipe V Freitas

\e[4mUsage\e[0m: 
reconstruct-species-tree.sh   --iqtree-out <path/to/iqtree-out> \\
                              --astral-in <path/to/astral-in> \\
                              --astral-out <path/to/astral-out> \\
                              --astral-jar <Astral.X.X.X.jar>

\e[4mRequired arguments\e[0m:
 --iqtree-out           Path to gene trees (produced by estimate-trees.sh)
 --astral-in            Path to store input of astral
 --astral-out           Path to store output of astral
 --astral-jar           Path to astral .jar file

\e[4mOptional arguments\e[0m:
 --only-by-gene         Only reconstruct species tree for 'concatenated-by-gene/' files
 --only-by-region       Only reconstruct species tree for 'concatenated-by-region/' files
 
 "
exit 2
}

# Option strings for arg parser
SHORT=h
LONG=help,astral-in:,iqtree-out:,astral-out:,astral-jar:,only-by-gene,only-by-region

# Set deafult values
ONLY_BY_GENE="False"
ONLY_BY_REGION="False"

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
		--astral-in )
		ASTRAL_IN="$2"
		shift 2
		;;
		--iqtree-out )
		IQTREE_OUT="$2"
		shift 2
		;;
		--astral-out )
		ASTRAL_OUT="$2"
		shift 2
		;;
		--astral-jar )
		ASTRAL_JAR="$2"
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

# Checking if any required args is empty
if [ -z "${ASTRAL_OUT}" ] || [ -z "${IQTREE_OUT}" ] || \
   [ -z "${ASTRAL_IN}" ] || [ -z "${ASTRAL_JAR}" ]; then
	error_exit "Please, supply all arguments correctly."
fi

BY_GENE="$ASTRAL_IN"/by-gene
BY_REGION="$ASTRAL_IN"/by-region

# prepare dir by gene
if [ "$ONLY_BY_REGION" == "False" ]; then
	echo "- Preparing species-tree by gene"
	mkdir -p $BY_GENE
	# only genic regions
	cat "$IQTREE_OUT"/concatenated-by-gene/*treefile \
		> "$BY_GENE"/only-genic-regions.tre
	# all regions (add intergenic)
	cat "$IQTREE_OUT"/intergenic-regions/*treefile \
		"$BY_GENE"/only-genic-regions.tre \
		> "$BY_GENE"/all-regions.tre
fi

# prepare dirs by region
if [ "$ONLY_BY_GENE" == "False" ]; then
	echo "- Preparing species-tree by region"
	mkdir -p $BY_REGION
	# only exons
	cat "$IQTREE_OUT"/concatenated-by-region/*__exon-*treefile \
		> "$BY_REGION"/only-exons.tre
	# only introns
	cat "$IQTREE_OUT"/concatenated-by-region/*__introns.treefile \
		> "$BY_REGION"/only-introns.tre
	# only genic regions (exons + introns)
	cat "$BY_REGION"/only-introns.tre "$BY_REGION"/only-exons.tre \
		> "$BY_REGION"/only-genic-regions.tre
	# all regions (exons + introns + intergenic regions)
	cat "$IQTREE_OUT"/intergenic-regions/*treefile \
		"$BY_REGION"/only-introns.tre "$BY_REGION"/only-exons.tre \
		> "$BY_REGION"/all-regions.tre
fi

echo "- Running ASTRAL"
# run by gene
if [ "$ONLY_BY_REGION" == "False" ]; then
	mkdir -p "$ASTRAL_OUT"/by-gene 
	echo "--- By gene, all regions"
	java -jar "$ASTRAL_JAR" \
		-i "$ASTRAL_IN"/by-gene/all-regions.tre \
		-o "$ASTRAL_OUT"/by-gene/all-regions.tre
	echo "--- By gene, genic regions"
	java -jar "$ASTRAL_JAR" \
		-i "$ASTRAL_IN"/by-gene/only-genic-regions.tre \
		-o "$ASTRAL_OUT"/by-gene/only-genic-regions.tre
fi

# run by region
if [ "$ONLY_BY_GENE" == "False" ]; then
	mkdir -p "$ASTRAL_OUT"/by-region
	echo "--- By region, all regions"
	java -jar "$ASTRAL_JAR" \
		-i "$ASTRAL_IN"/by-region/all-regions.tre \
		-o "$ASTRAL_OUT"/by-region/all-regions.tre 
	echo "--- By region, genic regions"
	java -jar "$ASTRAL_JAR" \
		-i "$ASTRAL_IN"/by-region/only-genic-regions.tre \
		-o "$ASTRAL_OUT"/by-region/only-genic-regions.tre 
	echo "--- By region, introns"
	java -jar "$ASTRAL_JAR" \
		-i "$ASTRAL_IN"/by-region/only-introns.tre \
		-o "$ASTRAL_OUT"/by-region/only-introns.tre
	echo "--- By region, exons"
	java -jar "$ASTRAL_JAR" \
		-i "$ASTRAL_IN"/by-region/only-exons.tre \
		-o "$ASTRAL_OUT"/by-region/only-exons.tre
fi

echo "- All done. Bye"