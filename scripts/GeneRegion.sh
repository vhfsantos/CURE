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
# "Bye" message function (green text)
BYEmsg() {
	echo -e "\033[1;32m[ CURE v$VERSION | $(date +%Y-%m-%d" "%H:%M:%S) ] ALL DONE!\033[0m"
	echo -e "
\033[1;32mThis is a goodbye message.
Probably will have citation, mail contact, and more.
Bye\033[0m"
}

# Usage function
Usage() {
echo -e "
-----------------------------------------------------------
CURE: an automated and parallel pipeline for UCE curation
-----------------------------------------------------------
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
  -f, --filter-string     UCEs whose name beggins with this string will be discarted (Default: "_")
  --only-by-gene          Concatenate UCEs only by gene. Concatenation by region will be discarted
                          with this flag
  --only-by-region        Concatenate UCEs only by region. Concatenation by gene will be discarted
                          with this flag"
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
	for app in $CONDA_PREFIX/bin/blastn $CONDA_PREFIX/bin/parallel $CONDA_PREFIX/bin/phyluce_align_concatenate_alignments; do
        	command -v $app >/dev/null 2>&1 || \
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
ONLY_BY_GENE="False"
ONLY_BY_REGION="False"
FILTER="_"

# Option strings for arg parser
SHORT=hb:r:g:p:o:f:t:
LONG=help,baits:,version:,reference:,gff:,phyluce-nexus:,output:,filter-string:,threads:,only-by-gene,only-by-region


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

# Checking if any required args is empty
if [ -z "${NEXUS_DIR}" ] || [ -z "${BAITS_FILE}" ] || \
	[ -z "${REFERENCE_GENOME}" ] || [ -z "${GFF}" ] || \
	[ -z "${OUTPUT}" ]; then
	error_exit "Please, supply all arguments correctly."
fi

# Checking if both flags were signed
if [ "$ONLY_BY_GENE" == "True" ] && [ "$ONLY_BY_REGION" == "True" ]; then
	error_exit "Please, chose one: --only-by-gene or --only-by-region"
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
ONLY_BY_REGION? $ONLY_BY_REGION
FILTER_STRING $FILTER
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

GENOME=$(basename "${REFERENCE_GENOME}")
UCE_KIT_SUMMARY="${OUTPUT}/uce_kit_output/${GENOME}.uce_kit_summary"

if [[ ! -f "${UCE_KIT_SUMMARY}" ]]; then
	log "Running uce_kit.py..."
	$CONDA_PREFIX/bin/python "${HOME_DIR}"/uce_kit/uce_kit.py run_pipeline \
		"${BAITS_FILE}" \
		"${REFERENCE_GENOME}" \
		"${GFF}" "${OUTPUT}"/uce_kit_output \
		-lines -merge -filt "$FILTER" > "${LOGDIR}"/uce_kit.log 2>&1
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
if [ ! -f "${OUTPUT}/CURE-exons.txt" ]; then
	log "Assinging UCEs to exons..."
	EXONS=$(mktemp)
	# Create temp files to count each assignment type
	ASSIGNED2EXON=${OUTPUT}/tmp/assign2ex
	ASSIGNED2INTRON=${OUTPUT}/tmp/assign2in
	ASSIGNED2EXON_NOTMISSING=${OUTPUT}/tmp/assign2exNotMissing
	ASSIGNED2INTRON_NOTMISSING=${OUTPUT}/tmp/assign2inNotMissing
	touch $ASSIGNED2EXON
	touch $ASSIGNED2INTRON
	BEFOREALL=${OUTPUT}/tmp/2Bassign
	# Get names of all UCEs available.
	find ${NEXUSCOPYex} -type f | xargs -L1 -I{} basename "{}" > "$BEFOREALL"
	TOTAL_UCEs=$(wc -l < $BEFOREALL)
	# Get only exons from uce_kit summary.
	awk \
		-v FS="\t" \
		-v OFS="\t" \
		'{ if( $4 == "exon" ){print $1, $2, $5, $6, $10} }' \
		"${UCE_KIT_SUMMARY}" | uniq \
		> "$EXONS"
	# Run ParseResults for exons.
	$CONDA_PREFIX/bin/python "${HOME_DIR}"/parse-results.py "$EXONS" \
	                "$FILTER" "${OUTPUT}"/CURE-exons.txt exon
	# Read output from script.
	# Lines are uce name, exon ID and gene ID.
	# So, copy uce name from nexus dir to gene/exon dir.
	while IFS=$'\t' read uce exonID geneID; do
		UCEfile=${NEXUSCOPYex}/${uce/-/_}
		GENEdir=${NEXUSTOCONCAT}/exons/$geneID
		EXONdir=${GENEdir}/$exonID
		echo "$uce" >> $ASSIGNED2EXON
		if [[ ! -f ${UCEfile} ]]; then
			warn "${uce%.*} already assigned or missing in nexus dir. \
			Skipping this file..."
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
	$CONDA_PREFIX/bin/python "${HOME_DIR}"/parse-results.py \
	        "$INTRONS" "$FILTER" "${OUTPUT}"/CURE-introns.txt intron
	# Same as previous 'while', but for introns.
	while IFS=$'\t' read uce geneID; do
		UCEfile=${NEXUSCOPYex}/${uce/-/_}
		GENEdir=${NEXUSTOCONCAT}/introns/$geneID
		echo "$uce" >> $ASSIGNED2INTRON
		if [[ ! -f "${UCEfile}" ]]; then
			warn "${UCEfile} already assigned or missing in nexus dir. \
			Skipping this file..."
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
	$CONDA_PREFIX/bin/python "${HOME_DIR}"/parse-results.py \
	        "$INTERGENIC" "$FILTER" "${OUTPUT}"/CURE-intergenic.txt \
	        intergenic
	# Same as previous 'while', but only with UCE name.
	while IFS=$'\t' read uce; do
		UCEfile=${NEXUSCOPYex}/${uce/-/_}
		if [[ ! -f ${UCEfile} ]]; then
			warn "${uce%.*} already assigned or missing in nexus dir. \
			Skipping this file..."
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
echo "CHEGOU ATÈ AQUIIIIIIIII"
grep -o -f $ASSIGNED2INTRON $BEFOREALL > $ASSIGNED2INTRON_NOTMISSING
grep -o -f $ASSIGNED2EXON $BEFOREALL > $ASSIGNED2EXON_NOTMISSING
# Count each type of assignment:
# E: exon
# I: intron
# EI: exon-and-intron
# N: intergenic regions
EI=$(grep -f $ASSIGNED2INTRON_NOTMISSING $ASSIGNED2EXON_NOTMISSING \
	| uniq | wc -l)
E=$(grep -v -f $ASSIGNED2INTRON_NOTMISSING $ASSIGNED2EXON_NOTMISSING \
	| uniq | wc -l)
I=$(grep -v -f $ASSIGNED2EXON_NOTMISSING $ASSIGNED2INTRON_NOTMISSING \
	| uniq | wc -l)
N=$(find $INTERGENIC_DIR -type f | wc -l)
# move unassigned UCEs to 'intergenic' dir
UN=$(ls ${NEXUSCOPYex} | wc -l)
echo "CHGEOU AQUI TBM "
mv ${NEXUSCOPYex}/* "$INTERGENIC_DIR"
echo "E AQUI TBM POOOOOO"
# print stats
log "--------------------------------------------"
log "----------------- SUMMARY ------------------"
log "--------------------------------------------"
log "$TOTAL_UCEs UCEs in NEXUS dir:"
log "$E assigned to exons"
log "$I assigned to introns"
log "$EI assigned to both (accounted for exons)"
log "$N assigned to intergenic regions"
log "$UN unassigned (accounted for intergenic)"
log "--------------------------------------------"
# Printing to csv
echo "type,uce.count
total,$TOTAL_UCEs
exons,$E
introns,$I
exon-and-intron,$EI
intergenic,$N
unassigned, $UN" > ${OUTPUT}/CURE-summary.csv

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
	log "Concatenating UCEs by region (exons)..."
	DIRS=$(find "${NEXUSTOCONCAT}"/exons/ -mindepth 2 -type d)
	N_DIRS=$(find "${NEXUSTOCONCAT}"/exons/ -mindepth 2 -type d | wc -l)
	AUX=0
	for dir in $DIRS; do
		# Output dir is named as gene-[geneID]__exon-[exonID]/
		# (separated by two underscores)
		GENEID__EXONID=$(echo "$dir" | rev \
			| cut -d "/" -f 1,2 | rev | sed 's/\//__/g' | \
			sed 's/-/_/g')
		GENEID=$(echo "$dir" | rev | cut -d "/" -f 2 | rev)
		AUX=$(( AUX + 1 ))
		# concatenate with phyluce in parallel
		$CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS" \
			$CONDA_PREFIX/bin/phyluce_align_concatenate_alignments \
			--alignments "$dir" --nexus --log "${LOGDIR}" \
			--output "${ALL_EXONS_DIR}"/"${GENEID__EXONID}" > /dev/null 2>&1
		"${HOME_DIR}"/progress-bar.sh $AUX "$N_DIRS"
	done
	$CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
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
	log "Concatenating UCEs by region (introns)..."
	DIRS=$(find "${NEXUSTOCONCAT}"/introns/ -mindepth 1 -type d)
	N_DIRS=$(find "${NEXUSTOCONCAT}"/introns/ -mindepth 1 -type d | wc -l)
	AUX=0
	for dir in $DIRS; do
		# the output dir will be named as gene-[geneID]
		GENEID=$(echo "$dir" | rev | cut -d "/" -f 1 | rev  | \
			sed 's/-/_/g')
		AUX=$(( AUX + 1 ))
		$CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS" \
			$CONDA_PREFIX/bin/phyluce_align_concatenate_alignments \
			--alignments "$dir" --nexus --log "${LOGDIR}" \
			--output "${ALL_INTRONS_DIR}"/"${GENEID}" > /dev/null 2>&1
		"${HOME_DIR}"/scripts/progress-bar.sh $AUX "$N_DIRS"
	done
	$CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
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
	PREPDIRS=$(find "${GENE_MERGED_PREP_BASEDIR}" -mindepth 1 -type d)
	N_PREPDIRS=$(find "${GENE_MERGED_PREP_BASEDIR}" -mindepth 1 -type d | wc -l)
	AUX=0
	for dir in $PREPDIRS; do
		GENEID=$(echo "$dir" | rev | cut -d "/" -f 1 | rev | \
			sed 's/-/_/g')
		AUX=$(( AUX + 1 ))
		$CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS" \
			$CONDA_PREFIX/bin/phyluce_align_concatenate_alignments \
			--alignments "$dir" \
			--phylip --input-format nexus --log "${LOGDIR}" \
			--output "${GENE_MERGED_BASEDIR}"/"${GENEID}" \
			> /dev/null 2> /dev/null
		"${HOME_DIR}"/scripts/progress-bar.sh $AUX "$N_PREPDIRS"
	done
	$CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
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
BY_REGION_DIR="$OUTPUT"/concatenated-by-region

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
	log "Writing alignments concatenated by region..."
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
rm -rf ${OUTPUT}/tmp
rm -rf ${OUTPUT}/phyluce_align_concatenate_alignments.log
DONEmsg
