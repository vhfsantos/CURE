#!/bin/bash
set -e
VERSION="$2"

trap ctrl_c SIGINT
trap 'echo - Ignoring a SIGHUP received..' SIGHUP

function ctrl_c() {
	warn "Canceling all parallel jobs..."
	killall parallel
	warn "Bye"
	exit 1
}
# Log message function (normal text)
log() {
	echo "[ CURE v$VERSION | $(date +%Y-%m-%d" "%H:%M:%S) ] $1"
}
# Warning message function (yellow text)
warn() {
	echo -e "\033[1;33m[ CURE v$VERSION | $(date +%Y-%m-%d" "%H:%M:%S) ] $1\033[0m"
	}
# "Done" message function
DONEmsg() {
	echo -e "[ CURE v$VERSION | $(date +%Y-%m-%d" "%H:%M:%S) ] Done!"
}

# Usage function
Usage() {
echo -e "
-----------------------------------------
CURE: \e[4mC\e[0muration of \e[4mU\e[0mltraconse\e[4mR\e[0mved \e[4mE\e[0mlements
-----------------------------------------
by Vinícius Franceshini-Santos & Felipe Freitas
version "$VERSION"

\e[4mUsage\e[0m:
 CURE GeneRegion --baits <baits.fa> --reference <reference.fasta> \\
                 --gff <reference.gff> --phyluce-nexus <nexus_dir> \\
                 --output <output_dir>

\e[4mRequired arguments\e[0m:
  -b, --baits             Path to UCE baits file

  -r, --reference         Path to reference genome file

  -g, --gff               Path to annotation file in gff format

  -p, --phyluce-nexus     Path to the directory with nexus files created with
                          phyluce

  -o, --output            Output directory

  \e[4mOptional arguments\e[0m:
  -t, --threads           Number of threads for the analysis (Default: 2)

  --min-identity          Minimum bait/genome identity, 0-100 (Default: 99)

  --min-length            Minimum aligned probe length (Default: 120)

  -f, --filter-string     UCEs whose name beggins with this string will be discarted (Default: "_")

  --only-by-gene          Concatenate UCEs only by gene. Concatenation by genic region will not
                          be performed with this flag
						  
  --only-by-genic-region  Concatenate UCEs only by genic region. Concatenation by gene will not 
                          be performed with this flag"
exit 2
}

# Error function
error_exit() {
	msg=$1
	>&2 echo -e "\033[1;31mERROR: ${msg}\033[0m"
	Usage
}

# Check dependencies
check_deps() {
	for app in blastn parallel phyluce_align_concatenate_alignments; do
        	command -v $CONDA_PREFIX/bin/$app >/dev/null 2>&1 || \
			error_exit "Cannot find ${app} in your PATH variable\nDid you activate the cure environment?"
	done
}

# Check deps
check_deps

# Print usage if no arg
if [ $# -eq 0 ]; then usage; fi

# Setting home dir for utils calling
tmp=$(realpath "$0")
HOME_DIR=${tmp%/*}
THREADS=2
MIN_IDENTITY=99
MIN_LENGTH=120
ONLY_BY_GENE="False"
ONLY_BY_REGION="False"
FILTER="_"

# Option strings for arg parser
SHORT=hb:r:g:p:o:f:t:
LONG=help,baits:,version:,reference:,gff:,phyluce-nexus:,output:,filter-string:,threads:,min-identity:,min-length:,only-by-gene,only-by-genic-region


# Read options
OPTS=$(getopt --options $SHORT --long $LONG --name "$0" -- "$@")
if [ $? != 0 ] ; then error_exit "Failed to parse options...exiting."; fi
eval set -- "$OPTS"

# Extracting arguments
while true ; do
	case "$1" in
		-h | --help )
		Usage
		;;
		--version )
		VERSION="$2"
		shift 2
		;;
		-p | --phyluce-nexus )
		NEXUS_DIR="$2"
		shift 2
		;;
		-b | --baits )
		BAITS_FILE="$2"
		shift 2
		;;
		-r | --reference )
		REFERENCE_GENOME="$2"
		shift 2
		;;
		-t | --threads )
		THREADS="$2"
		shift 2
		;;
		--min-identity )
		MIN_IDENTITY="$2"
		shift 2
		;;
		--min-length )
		MIN_LENGTH="$2"
		shift 2
		;;
		-g | --gff )
		GFF="$2"
		shift 2
		;;
		-o | --output )
		OUTPUT="$2"
		shift 2
		;;
		-f | --filter-string )
		FILTER="$2"
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
		-- )
		shift
		break
		;;
		*)
		error_exit "Please, supply all arguments correctly."
		;;
	esac
done

# Checking if any required args is empty
if [ -z "${NEXUS_DIR}" ] || [ -z "${BAITS_FILE}" ] || \
	[ -z "${REFERENCE_GENOME}" ] || [ -z "${GFF}" ] || \
	[ -z "${OUTPUT}" ]; then
	error_exit "Please, supply all arguments correctly."
fi

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || error_exit "--threads must be a positive integer"
[[ "$MIN_LENGTH" =~ ^[1-9][0-9]*$ ]] || error_exit "--min-length must be a positive integer"
[[ "$MIN_IDENTITY" =~ ^([0-9]+)(\.[0-9]+)?$ ]] || error_exit "--min-identity must be between 0 and 100"
awk -v n="$MIN_IDENTITY" 'BEGIN { exit !(n >= 0 && n <= 100) }' || error_exit "--min-identity must be between 0 and 100"
GENOME=$(basename "$REFERENCE_GENOME")
UCE_KIT_SUMMARY="${OUTPUT}/uce_kit_output/${GENOME}.uce_kit_summary"
UCE_KIT_PARAMETERS="${OUTPUT}/uce_kit_output/parameters.txt"
EXPECTED_PARAMETERS=$(printf '%s\n' "$BAITS_FILE" "$REFERENCE_GENOME" "$GFF" "$NEXUS_DIR" \
	"$FILTER" "$MIN_IDENTITY" "$MIN_LENGTH" "$ONLY_BY_GENE" "$ONLY_BY_REGION")
if [[ -f "$UCE_KIT_SUMMARY" ]]; then
	[[ -f "$UCE_KIT_PARAMETERS" ]] || error_exit "Existing uce_kit results have no parameter record; use a new output directory"
	[[ "$(cat "$UCE_KIT_PARAMETERS")" == "$EXPECTED_PARAMETERS" ]] || error_exit "Existing uce_kit results used different inputs or thresholds; use a new output directory"
fi
if [[ -f "$OUTPUT/GeneRegion.complete" ]]; then
	log "GeneRegion already completed in $OUTPUT"
	exit 0
fi
if [[ -d "$OUTPUT/tmp" || -f "$OUTPUT/CURE_stats.csv" ]]; then
	warn "Previous GeneRegion run was incomplete; restarting classification"
	rm -rf -- "$OUTPUT/tmp" "$OUTPUT/concatenated-by-gene" \
		"$OUTPUT/concatenated-by-genic-region" "$OUTPUT/intergenic-regions" \
		"$OUTPUT/unassigned-regions" "$OUTPUT/unclassified-mapped"
	rm -f -- "$OUTPUT/CURE-exons.txt" "$OUTPUT/CURE-introns.txt" \
		"$OUTPUT/CURE-intergenic.txt" "$OUTPUT/CURE_stats.csv" "$OUTPUT/CURE_stats.pdf"
fi

# Checking if both flags were signed
if [ "$ONLY_BY_GENE" == "True" ] && [ "$ONLY_BY_REGION" == "True" ]; then
	error_exit "Please, chose one: --only-by-gene or --only-by-genic-region"
fi

# Print parameters for debuggin
echo "=========================================================================
                             CURE (v$VERSION)
          An automated and parallel pipeline for UCE curation
          by Vinícius Franceschini-Santos & Felipe V Freitas
-------------------------------------------------------------------------
                              GeneRegion
-------------------------------------------------------------------------
NEXUS: ${NEXUS_DIR}
BAITS: ${BAITS_FILE}
GENOME: ${REFERENCE_GENOME}
GFF: ${GFF}
OUTDIR: ${OUTPUT}
THREADS: $THREADS
ONLY_BY_GENE? $ONLY_BY_GENE
ONLY_BY_GENIC_REGION? $ONLY_BY_REGION
FILTER_STRING $FILTER
MIN_IDENTITY $MIN_IDENTITY
MIN_LENGTH $MIN_LENGTH
-------------------------------------------------------------------------"

#=============================================================
#====                       STEP 0:                       ====
#====                 Preparing input data                ====
#=============================================================

# Create a copy os nexus dir.
# Allow me to 'move' UCEs from there.
# Moving avoids one UCE assigned to multiple exons/introns.
# Also, if UCE is placed on one exon and one intron, it
# is assigned to the exon.
NEXUSCOPYex=${OUTPUT}/tmp/000-nexus-copy/
mkdir -p ${NEXUSCOPYex}

if [ -z "$(ls -A "${NEXUSCOPYex}")" ]; then
	log "Preparing input data..."
	$CONDA_PREFIX/bin/parallel \
	--will-cite --max-procs ${THREADS} \
	cp {} ${NEXUSCOPYex} ::: ${NEXUS_DIR}/* 2> /dev/null
else
	warn "Input preparation already done. Skipping"
fi

#=============================================================
#====                       STEP 1:                       ====
#====         Run a customize version of uce_kit.py       ====
#=============================================================
# Original:
# https://github.com/calacademy-research/ccgutils/tree/master/uce_types

mkdir -p "${OUTPUT}"/tmp/
mkdir -p "${OUTPUT}"/uce_kit_output
LOGDIR="${OUTPUT}"/logfiles/
mkdir -p "${LOGDIR}"

if [[ ! -f "${UCE_KIT_SUMMARY}" ]]; then
	printf '%s\n' "$EXPECTED_PARAMETERS" > "$UCE_KIT_PARAMETERS"
	log "Running uce_kit.py..."
	$CONDA_PREFIX/bin/python "${HOME_DIR}"/uce_kit/uce_kit.py run_pipeline \
		"${BAITS_FILE}" \
		"${REFERENCE_GENOME}" \
		"${GFF}" "${OUTPUT}"/uce_kit_output \
		-lines -merge -filt "$FILTER" -pct "$MIN_IDENTITY" -len "$MIN_LENGTH" > "${LOGDIR}"/uce_kit.log 2>&1
	[[ -s "$UCE_KIT_SUMMARY" ]] || error_exit "uce_kit produced no classification summary; see ${LOGDIR}/uce_kit.log"
	DONEmsg
	else
		warn "Script uce_kit.py already run. Skipping..."
fi

#=============================================================
#====                       STEP 2:                       ====
#====  Organize UCEs according to their characterization  ====
#=============================================================

# I treat each type differently.
# Let's beggin with exons.
EXONCOUNT=0
NEXUSTOCONCAT=${OUTPUT}/tmp/00-characterized-uces/
mkdir -p "${NEXUSTOCONCAT}/exons" "${NEXUSTOCONCAT}/introns"
if [ ! -f "${OUTPUT}/CURE-exons.txt" ]; then
	log "Assinging UCEs to exons..."
	EXONS=$(mktemp)
	# Create temp files to count each assignment type
	#ASSIGNED2EXON=${OUTPUT}/tmp/assign2ex
	#ASSIGNED2INTRON=${OUTPUT}/tmp/assign2in
	#ASSIGNED2EXON_NOTMISSING=${OUTPUT}/tmp/assign2exNotMissing
	#ASSIGNED2INTRON_NOTMISSING=${OUTPUT}/tmp/assign2inNotMissing
	#touch $ASSIGNED2EXON
	#touch $ASSIGNED2INTRON
	#BEFOREALL=${OUTPUT}/tmp/2Bassign
	# Get names of all UCEs available.
	#find ${NEXUSCOPYex} -type f | xargs -L1 -I{} basename "{}" > "$BEFOREALL"
	#TOTAL_UCEs=$(wc -l < $BEFOREALL)
	# Get only exons from uce_kit summary.
	awk \
		-v FS="\t" \
		-v OFS="\t" \
		'{ if( $4 == "exon" ){print $1, $2, $5, $6, $10} }' \
		"${UCE_KIT_SUMMARY}" | uniq \
		> "$EXONS"
	# Run ParseResults for exons.
	$CONDA_PREFIX/bin/python "${HOME_DIR}"/parse_results.py "$EXONS" \
	                "$FILTER" "${OUTPUT}"/CURE-exons.txt exon "${GFF}"
	# Read output from script.
	# Lines are uce name, exon ID and gene ID.
	# So, copy uce name from nexus dir to gene/exon dir.
	while IFS=$'\t' read uce exonID geneID; do
		UCEfile=${NEXUSCOPYex}/${uce}
		GENEdir=${NEXUSTOCONCAT}/exons/$geneID
		EXONdir=${GENEdir}/$exonID
		if [[ ! -f ${UCEfile} ]]; then
			#warn "${uce%.*} already assigned or missing in nexus dir. \
			#Skipping this file..."
			continue
		else
			mkdir -p "$EXONdir"
			EXONCOUNT=$(( $EXONCOUNT+1 ))
			mv "${UCEfile}" "${EXONdir}"/"$uce"
		fi

	done < "${OUTPUT}"/CURE-exons.txt
	DONEmsg
else
	warn "Already assigned UCEs to exons. Skipping..."
fi

# Doing (almost) the same for introns.
# Introns do not have IDs like exons.
# Output dir will be only geneID.
if [ ! -f "${OUTPUT}/CURE-introns.txt" ]; then
	log "Assigning UCEs to introns..."
	INTRONS=$(mktemp)
	awk \
		-v FS="\t" \
		-v OFS="\t" \
		'{ if( $4 == "intron" ){print $1, $2, $5, $6, $10} }' \
		"${UCE_KIT_SUMMARY}" | uniq \
		> "$INTRONS"
	# Run ParseResults for introns.
	$CONDA_PREFIX/bin/python "${HOME_DIR}"/parse_results.py \
	        "$INTRONS" "$FILTER" "${OUTPUT}"/CURE-introns.txt intron
	# Same as previous 'while', but for introns.
	while IFS=$'\t' read uce geneID; do
		UCEfile=${NEXUSCOPYex}/${uce}
		GENEdir=${NEXUSTOCONCAT}/introns/$geneID
#		echo "$uce" >> $ASSIGNED2INTRON
		if [[ ! -f "${UCEfile}" ]]; then
			#warn "${UCEfile} already assigned or missing in nexus dir. \
			#Skipping this file..."
			continue
		else
			mkdir -p "$GENEdir"
			mv "${UCEfile}" "${GENEdir}"/"$uce"
		fi
	done < "${OUTPUT}"/CURE-introns.txt
	DONEmsg
else
	warn "Already assigned UCEs to introns. Skipping..."
fi

# Now, the same but simpler for intergenic regions.
# Simpler because I don't need geneID here.
if [ ! -f "${OUTPUT}/CURE-intergenic.txt" ]; then
	log "Assigning UCEs to intergenic regions..."
	INTERGENIC_DIR=${OUTPUT}/intergenic-regions
	mkdir -p "$INTERGENIC_DIR"
	INTERGENIC=$(mktemp)
	awk \
		-v FS="\t" \
		-v OFS="\t" \
		'{ if( $6 == "intergenic" ){print $1} }' \
		"${UCE_KIT_SUMMARY}" | uniq \
		> "$INTERGENIC"
	# Run ParseResults for intergenic region.
	$CONDA_PREFIX/bin/python "${HOME_DIR}"/parse_results.py \
	        "$INTERGENIC" "$FILTER" "${OUTPUT}"/CURE-intergenic.txt \
	        intergenic
	# Same as previous 'while', but only with UCE name.
	while IFS=$'\t' read uce; do
		UCEfile=${NEXUSCOPYex}/${uce}
		if [[ ! -f ${UCEfile} ]]; then
			#warn "${uce%.*} already assigned or missing in nexus dir. \
			#Skipping this file..."
			continue
		else
			mv "${UCEfile}" "${INTERGENIC_DIR}"/"$uce"
		fi
	done < "${OUTPUT}"/CURE-intergenic.txt
	DONEmsg
else
	warn "Already assigned UCEs to intergenic regions. Skipping..."
fi
# Compute statistics.
# I didn't use uce_kit stats because its based on UCEs in baits file.
# A lot of UCEs are in baits file but not in NEXUS dir.
# So real stats for this should be based on UCEs in NEXUS dir.
UCE_IN_INTER=$(mktemp)
UCE_IN_EXONS=$(mktemp)
UCE_IN_INTRONS=$(mktemp)
ALL_UCES=$(mktemp)

cat "${OUTPUT}/CURE-intergenic.txt" >> $UCE_IN_INTER
cat "${OUTPUT}/CURE-exons.txt" | cut -f1 | uniq >> $UCE_IN_EXONS
cat "${OUTPUT}/CURE-introns.txt" | cut -f1 | uniq >> $UCE_IN_INTRONS
find "$NEXUS_DIR" -maxdepth 1 -type f -name '*.nexus' -printf '%f\n' > "$ALL_UCES"

# plot venn diagram and get statistics

"${CONDA_PREFIX}/bin/python" "${HOME_DIR}/CURE_GeneRegion_plot_venn.py" \
	"$UCE_IN_INTER" "$UCE_IN_EXONS" "$UCE_IN_INTRONS" "$ALL_UCES" "$OUTPUT" "$VERSION" \
	"${OUTPUT}/uce_kit_output/$(basename "${REFERENCE_GENOME%.*}")_uce_locations.tsv"
"${CONDA_PREFIX}/bin/python" "${HOME_DIR}/gene_region_remainder.py" \
	"$NEXUSCOPYex" "$OUTPUT" \
	"${OUTPUT}/uce_kit_output/$(basename "${REFERENCE_GENOME%.*}")_uce_locations.tsv"

#=============================================================
#====                       STEP 3:                       ====
#====        Concatenating UCEs by region (exons)         ====
#=============================================================

ALL_EXONS_DIR=${OUTPUT}/tmp/01-exon-data/
mkdir -p "$ALL_EXONS_DIR"

# This step turn this:
# - output/tmp/00-.../
#    - gene-X/
#        - exon-1/
#        	 - uce-1.nexus
#            - uce-2.nexus
#
# Into this:
# - output/tmp/01-.../exons/
#    - gene-X__exon-1/
#        - gene-X__exon-1.nexus

if [ -z "$(ls -A "${ALL_EXONS_DIR}")" ]; then
	log "Merging UCEs in the same exons..."
	"$CONDA_PREFIX/bin/python" "${HOME_DIR}/gene_region_concat.py" exon "${NEXUSTOCONCAT}/exons" "$ALL_EXONS_DIR" "${LOGDIR}/exon-concatenation" --threads "$THREADS"
	DONEmsg
else
	warn "Already concatenated UCEs from same exons. Skipping..."
fi

#=============================================================
#====                       STEP 4:                       ====
#====        Concatenating UCEs by region (introns)       ====
#=============================================================

ALL_INTRONS_DIR=${OUTPUT}/tmp/02-intron-data/
mkdir -p "$ALL_INTRONS_DIR"

# This step turn this:
# - output/tmp/00-.../introns/
#    - gene-X/
#       - uce-1.nexus
#       - uce-2.nexus
#
# Into this:
# - output/tmp/02-.../
#    - gene-X/
#        - gene-X.nexus

if [ -z "$(ls -A "${ALL_INTRONS_DIR}")" ]; then
	log "Merging UCEs in introns of the same gene..."
	"$CONDA_PREFIX/bin/python" "${HOME_DIR}/gene_region_concat.py" intron "${NEXUSTOCONCAT}/introns" "$ALL_INTRONS_DIR" "${LOGDIR}/intron-concatenation" --threads "$THREADS"
	DONEmsg
else
	warn "Already concatenate UCEs from same introns. Skipping..."
fi

#=============================================================
#====                       STEP 5:                       ====
#====             Prepare exons to cat by gene            ====
#=============================================================

GENE_MERGED_PREP_BASEDIR=${OUTPUT}/tmp/03-pre-concat/
mkdir -p "$GENE_MERGED_PREP_BASEDIR"

# This step turns this:
# - output/tmp/01-.../exons/
#    - gene-X__exon-1/
#        - gene-X__exon-1.nexus
#    - gene-X__exon-2/
#        - gene-X__exon-2.nexus
#
# Into this:
# - output/tmp/03-.../
#    - gene-X/
#        - exon-1.nexus
#        - exon-2.nexus

if [ -z "$(ls -A "${GENE_MERGED_PREP_BASEDIR}")" ]; then
	DIRS=$(find "${NEXUSTOCONCAT}"/exons/ -mindepth 2 -type d)
	N_DIRS=$(find "${NEXUSTOCONCAT}"/exons/ -mindepth 2 -type d | wc -l)
	log "Preparing exons for concatenation by gene..."
	for dir in $DIRS; do
		# Create a copy of it in another folder.
		# In this new location, the folders will be geneID.
		# Will contain all exononic UCEs (same gene) merged.
		EXONID=$(echo "$dir" | rev | cut -d "/" -f 1 | rev| \
			sed 's/-/_/g')
		GENEID=$(echo "$dir" | rev | cut -d "/" -f 2 | rev| \
			sed 's/-/_/g')
		GENEID_DIR=${GENE_MERGED_PREP_BASEDIR}/${GENEID}
		mkdir -p "${GENEID_DIR}"
		EXONIC_UCEs_PHYLIP=${ALL_EXONS_DIR}/${GENEID}__${EXONID}/${GENEID}__${EXONID}.nexus
		# Remove charsets from nexus to avoid error on next step.
		cat "$EXONIC_UCEs_PHYLIP" \
			| sed '/begin sets/,/end/d' \
			> "${GENEID_DIR}"/"${EXONID}".nexus
	done
	DONEmsg
else
	warn "Exons already prepared for concatenation by gene. Skipping..."
fi

#=============================================================
#====                       STEP 6:                       ====
#====            Prepare introns to cat by gene           ====
#=============================================================

# This step turns this:
# - output/tmp/03-.../
#    - gene-X/
#        - exon-1.nexus
#        - exon-2.nexus
# Into this:
# - output/tmp/03-.../
#    - gene-X/
#        - exon-1.nexus
#        - exon-2.nexus
#        - introns.nexus

if [ ! -f ${OUTPUT}/tmp/int.done ]; then
	DIRS=$(find "${NEXUSTOCONCAT}"/introns/ -mindepth 1 -type d)
	N_DIRS=$(find "${NEXUSTOCONCAT}"/introns/ -mindepth 1 -type d | wc -l)
	log "Preparing introns for concatenation by gene..."
	for dir in $DIRS; do
		GENEID=$(echo "$dir" | rev | cut -d "/" -f 1 | rev | \
			sed 's/-/_/g')
		GENEID_DIR=${GENE_MERGED_PREP_BASEDIR}/${GENEID}
		mkdir -p "${GENEID_DIR}"
		INTRONIC_UCEs_PHYLIP="${ALL_INTRONS_DIR}"/"${GENEID}"/"${GENEID}".nexus
		# remove the charsets from nexus file to avoid error on next step
		cat "$INTRONIC_UCEs_PHYLIP" \
			| sed '/begin sets/,/end/d' \
			> "${GENEID_DIR}"/introns.nexus
		touch ${OUTPUT}/tmp/int.done
	done
	DONEmsg
else
	warn "Introns already prepared for concatenation by gene. Skipping..."
fi

#=============================================================
#====                       STEP 7:                       ====
#====               Concatenating by gene                 ====
#=============================================================

GENE_MERGED_BASEDIR=${OUTPUT}/tmp/04-cat
mkdir -p "$GENE_MERGED_BASEDIR"

if [ -z "$(ls -A "${GENE_MERGED_BASEDIR}")" ]; then
	log "Concatenating UCEs by gene..."
	"$CONDA_PREFIX/bin/python" "${HOME_DIR}/gene_region_concat.py" gene "$GENE_MERGED_PREP_BASEDIR" "$GENE_MERGED_BASEDIR" "${LOGDIR}/gene-concatenation" --threads "$THREADS"
	DONEmsg
else
	warn "Concatenation by gene already done. Skipping..."
fi

#=============================================================
#====                       STEP 8:                       ====
#====             Removing remporary files                ====
#=============================================================

# Final output dirs
BY_GENE_DIR="$OUTPUT"/concatenated-by-gene
BY_REGION_DIR="$OUTPUT"/concatenated-by-genic-region

# cleaning tmp dirs
rm -rf ${NEXUSCOPYex} ${NEXUSTOCONCAT} "$ALL_EXONS_DIR" "$ALL_INTRONS_DIR"

# Create output concatenated by gene
if [ "$ONLY_BY_REGION" == "False" ]; then
	log "Writing alignments concatenated by gene..."
	mkdir -p "$BY_GENE_DIR"
	CATDIRS=$(find "$GENE_MERGED_BASEDIR" -mindepth 1 -type d)
	# If I kept final alignments in nexus format, charsets would
	# be merged inside the file by phyluce's script.
	# So I kept them in phylip format.
	for dir in $CATDIRS; do
		GENEID=$(echo "$dir" | rev | cut -d "/" -f 1 | rev)
		PHYLIP="$dir"/"$GENEID".phylip
		CHARSET="$dir"/"$GENEID".charsets
		NEWCHARSET="$BY_GENE_DIR"/"$GENEID".charsets

		# Moving alignment to final output.
		# Keep filename.
		mv "$PHYLIP" "$BY_GENE_DIR"

		# Moving chartset to final output.
		# Add #NEXUS to first line and remove "partition combined" line.
		sed '1i #NEXUS' "$CHARSET" | sed '/charpartition/d' > "$NEWCHARSET"
	done
fi

# Create output concatenated by region
if [ "$ONLY_BY_GENE" == "False" ]; then
	log "Writing alignments concatenated by genic region..."
	mkdir -p "$BY_REGION_DIR"
	PRECATFILES=$(find "$GENE_MERGED_PREP_BASEDIR" -type f)
	# Alignments format for 'by Gene' concatenation will be nexus.
	# A good way to avoid mixing alignments between approaches.
	for file in $PRECATFILES; do
		# Output nexus name will be geneID and regionID separated by
		# double underscore (__).
		# Exons ID already have single underscore.
		# (e.g. exon-XM_006561541.3-5)
		# So if I need to split geneID from regionID, I do it by "__".
		GENEID__REGIONID=$(echo $file | rev | cut -d "/" -f 1,2 \
			| rev | sed 's/\//__/g')
		NEWNEXUS="$BY_REGION_DIR"/"$GENEID__REGIONID"

		# Moving alignment to final output.
		mv "$file" "$NEWNEXUS"
	done
	DONEmsg
fi

log "Removing temporary files..."
for result_dir in "$BY_GENE_DIR" "$BY_REGION_DIR" "$INTERGENIC_DIR" \
	"$OUTPUT/unassigned-regions" "$OUTPUT/unclassified-mapped"; do
	if [[ -d "$result_dir" ]] && [[ -n "$(find "$result_dir" -type f -empty -print -quit)" ]]; then
		error_exit "An empty result was found in $result_dir; see $LOGDIR"
	fi
done
rm -rf -- "$OUTPUT/tmp"
rm -f -- "$OUTPUT/phyluce_align_concatenate_alignments.log"
touch "$OUTPUT/GeneRegion.complete"
DONEmsg
echo "- All done!

If you make use of the GeneRegion strategy in your research, please cite also the original paper describing the method:

Van Dam et al. (2021). Genomic Characterization and Curation of UCEs
Improves Species Tree Reconstruction. Systematic Biology 70(2):307-321.
https://doi.org/10.1093/sysbio/syaa063"
