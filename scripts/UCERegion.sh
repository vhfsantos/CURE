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

Usage() {
echo -e "
-----------------------------------------
CURE: \e[4mC\e[0muration of \e[4mU\e[0mltraconse\e[4mR\e[0mved \e[4mE\e[0mlements
-----------------------------------------
by Vinícius Franceshini-Santos & Felipe Freitas
version "$VERSION"

\e[4mUsage\e[0m:
 CURE UCERegion  --phyluce-nexus <nexus_dir> --output <output_dir>

\e[4mRequired arguments\e[0m:
  -p, --phyluce-nexus     Path to the directory with nexus files created with phyluce

  -o, --output            Output directory

\e[4mOptional arguments\e[0m:
  -t, --threads           Number of threads for the analysis (Default: 2)"

exit 2
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

# Error function
error_exit() {
       msg=$1
       >&2 echo -e "\033[1;31mERROR: ${msg}\033[0m"
       Usage
}

# Setting default variables
tmp=$(realpath "$0")
HOME_DIR=${tmp%/*}
THREADS=2
# Now SWSC is distributed by us!
SWSC_PATH="$HOME_DIR/SWSC_EN/SWSCEN.py"

# Option strings for arg parser
SHORT=hp:o:t:
LONG=help,phyluce-nexus:,output:,threads:,version:


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
        -t | --threads )
		THREADS="$2"
		shift 2
		;;
		-o | --output )
		OUTPUT="$2"
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

# Checking if any required args is empty
if [ -z "${NEXUS_DIR}" ] || [ -z "${OUTPUT}" ]; then
	error_exit "Please, supply all arguments correctly."
fi


# Print parameters for debuggin
echo "=========================================================================
                             CURE (v$VERSION)
          An automated and parallel pipeline for UCE curation
          by Vinícius Franceschini-Santos & Felipe V Freitas
-------------------------------------------------------------------------
                               UCERegion
-------------------------------------------------------------------------
NEXUS: ${NEXUS_DIR}
OUTDIR: ${OUTPUT}
THREADS: $THREADS
-------------------------------------------------------------------------"

#=============================================================
#====                       STEP 0:                       ====
#====                 Preparing input data                ====
#=============================================================

NEXUSCOPY=${OUTPUT}/tmp/000-nexus-copy/
mkdir -p ${NEXUSCOPY}
UCEs_in_subgroup=15
subgroups_dir="${OUTPUT}/tmp/001-subgroups"
N_UCES=$(find ${NEXUS_DIR} -maxdepth 1 -type f | wc -l)
n_subgroups=$((`find ${NEXUS_DIR} -maxdepth 1 -type f | wc -l`/16))

# workaround when less than 16 nexus:
if (( $n_subgroups == 0 )); then 
	n_subgroups=1
fi


if [ -z "$(ls -A "${NEXUSCOPY}")" ]; then
	log "Preparing input data..."
	log "Splitted $N_UCES UCEs into $n_subgroups subgroups"
        #  (1) Create a copy os nexus dir. Allows me to 'move' UCEs from there
	$CONDA_PREFIX/bin/parallel \
	--will-cite --max-procs ${THREADS} \
	cp {} ${NEXUSCOPY} ::: ${NEXUS_DIR}/* 2> /dev/null
        # (2) Split in subgroups of 15 UCEs for ease parallelization
        # thanks: https://stackoverflow.com/questions/29116212/split-a-folder-into-multiple-subfolders-in-terminal-bash-script

        for i in `seq 1 ${n_subgroups}`; do
                # create dir for subgroup
                mkdir -p "${subgroups_dir}/${i}"
                # copy files to there (replace uce-9.nexus > uce_9.nexus)
                find ${NEXUSCOPY} -type f | head -n $UCEs_in_subgroup \
                        | xargs -i mv "{}" "${subgroups_dir}/${i}"
        done
        DONEmsg
else
	warn "Input preparation already done. Skipping"
fi

#=============================================================
#====                       STEP 1:                       ====
#====               Concatenate UCEs for SWSC             ====
#=============================================================

SUBGROUPS_CAT=${OUTPUT}/tmp/002-subgroups-cat/
mkdir -p ${SUBGROUPS_CAT}

LOGDIR=${OUTPUT}/logfiles/
mkdir -p "${LOGDIR}"

if [ -z "$(ls -A "${SUBGROUPS_CAT}")" ]; then
	log "Concatenating UCEs subgroups for SWSC..."
        # concatenate with phyluce in parallel
        for sg in $(seq 1 $n_subgroups); do
                $CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS" \
                        $CONDA_PREFIX/bin/phyluce_align_concatenate_alignments \
                        --alignments "${subgroups_dir}/$sg" \
                        --nexus --log "${LOGDIR}" \
                        --output "${SUBGROUPS_CAT}/$sg" > /dev/null 2>&1
                bash "${HOME_DIR}"/progress-bar.sh $sg "$n_subgroups"
        done
        $CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
        DONEmsg
else
	warn "UCEs concatenation already done. Skipping"
fi


#=============================================================
#====                       STEP 2:                       ====
#====                      Run SWSC                       ====
#=============================================================

SWSC=${OUTPUT}/tmp/003-swsc/
mkdir -p ${SWSC}

if [ -z "$(ls -A "${SWSC}")" ]; then
	log "Running SWSC..."
        # run SWSC in parallel
        for sg in $(seq 1 $n_subgroups); do
                # (1) fix uce names in .nexus files
                sed -i 's/uce-/uce_/g' "${SUBGROUPS_CAT}/$sg/$sg.nexus"

                # (2) run SWSC
                $CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS" \
                        $CONDA_PREFIX/bin/python $SWSC_PATH \
                        $( realpath ${SUBGROUPS_CAT}/${sg}/${sg}.nexus ) \
                        $( realpath $SWSC ) > "${LOGDIR}"/swsc.log 2>&1
                bash "${HOME_DIR}"/progress-bar.sh $sg "$n_subgroups"
        done
        $CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
        DONEmsg
else
	warn "SWSC already run. Skipping"
fi

#=============================================================
#====                       STEP 3:                       ====
#====                   Parsing Results                   ====
#=============================================================

SWSC_PARSE=${OUTPUT}/tmp/004-swsc-parse
mkdir -p ${SWSC_PARSE}

# function to parse results
SWSCParser(){
	subgroup=$1
	SWSC=$2
	SUBGROUPS_CAT=$3
	SWSC_PARSE=$4
	OUTPUT=$5
	# (1) extract flanks coordenates for this subgroup
	   # 1st sed: adds 'charset' at the beggining of each line
	   # 2nd sed: adds 'begin sets;' as first line
	   # 3rd sed: adds 'end;' as last line
	grep '_right\|_core\|_left\|_all' ${SWSC}/${subgroup}.nexus_entropy_partition_finder.cfg \
		| sed 's/^/charset /g' \
		| sed '1 i\begin sets\;' \
		| sed -e '$aend;' > ${SWSC_PARSE}/${subgroup}.charsets

	# (2) remove old charsets from nexus file of this subgroup
	   # sed: delete all line after 'begin sets;'
	cat ${SUBGROUPS_CAT}/${subgroup}/${subgroup}.nexus \
		| sed '/begin sets;/,$d' > ${SWSC_PARSE}/${subgroup}.alignment

	# (3) create the new nexus file with flanks as charsets
	cat ${SWSC_PARSE}/${subgroup}.alignment \
		${SWSC_PARSE}/${subgroup}.charsets \
		> ${SWSC_PARSE}/${subgroup}.nexus

	# (4) call phyluce to split the nexus using the charsets
	$CONDA_PREFIX/bin/phyluce_align_split_concat_nexus_to_loci \
		--nexus ${SWSC_PARSE}/${subgroup}.nexus \
		--output ${SWSC_PARSE}/${subgroup}/ --log-path ${OUTPUT}/tmp/\
		--output-format nexus
}

# Need to export this function to call it with Parallel
export -f SWSCParser

if [ -z "$(ls -A "${SWSC_PARSE}")" ]; then
	log "Parsing results..."
        for sg in $(seq 1 $n_subgroups); do
		# calls function in parallel
                $CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS" \
                        SWSCParser $sg ${SWSC} ${SUBGROUPS_CAT} ${SWSC_PARSE} \
						${OUTPUT}> "${LOGDIR}"/swsc_parser.log 2>&1
                bash "${HOME_DIR}"/progress-bar.sh $sg "$n_subgroups"
        done
        $CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
        DONEmsg
else
	warn "SWSC results already parsed. Skipping"
fi

#=============================================================
#====                       STEP 4:                       ====
#====        Generating concat files and PF2 input        ====
#=============================================================

UCES_CAT=${OUTPUT}/tmp/005-concatenated-uces/
mkdir -p ${UCES_CAT}


if [ -z "$(ls -A "${UCES_CAT}")" ]; then
	log "Generating concatenated files and PF2 input..."
	# (1) Place all subgroups nexus to a same dir
	for sg in $(seq 1 $n_subgroups); do
	$CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS"
		cp ${SWSC_PARSE}/${sg}/* ${UCES_CAT}
	done
        $CONDA_PREFIX/bin/sem --will-cite --id $$ --wait

	# (2) Concatenate
	$CONDA_PREFIX/bin/phyluce_align_concatenate_alignments \
		--alignments "${UCES_CAT}" \
		--phylip --log "${LOGDIR}" \
		--output "${UCES_CAT}/concatenated-uces" > /dev/null 2>&1

	# (3) Header of .cfg file
	echo "## ALIGNMENT FILE ##
alignment = concatenated-uces-swscen.phylip;

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai <list> ##
models = GTR+G;

# MODEL SELECCTION: AIC | AICc | BIC #
model_selection = aicc;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]" > ${UCES_CAT}/concatenated-uces/header

	# (4) Body of .cfg file
	   # 1st sed: delete all after charpartition combined
	   # 2nd sed: remove 'begin sets;'
	   # 3rd sed: remove 'charset ' from beggining of lines
	cat "${UCES_CAT}"/concatenated-uces/concatenated-uces.charsets \
		| sed '/charpartition combined/,$d' \
		| sed 's/begin sets;//g' \
		| sed 's/^charset //g' | tail -n +3 \
		> ${UCES_CAT}/concatenated-uces/body

	# (5) Footer of .cfg file
	echo "
## SCHEMES, search: all | user | greedy | rcluster | hcluster | kmeans ##
[schemes]
search = rclusterf;" > ${UCES_CAT}/concatenated-uces/footer

	mkdir -p "${OUTPUT}/concatenated-uces/"
	cat ${UCES_CAT}/concatenated-uces/header \
		${UCES_CAT}/concatenated-uces/body \
		${UCES_CAT}/concatenated-uces/footer \
		> "${OUTPUT}/concatenated-uces/PF2-input.cfg"

	# (6) Get alignment from tmp dir
	mv ${UCES_CAT}/concatenated-uces/concatenated-uces.phylip "${OUTPUT}/concatenated-uces/concatenated-uces-swscen.phylip"



	### 1') Create charsets for IQ-tree
	cat "${UCES_CAT}"/concatenated-uces/concatenated-uces.charsets \
		| sed '1 i\#nexus' \
		| sed '/charpartition combined/d' > "${OUTPUT}/concatenated-uces/concatenated-uces-swscen.charsets"

else
	warn "PF2 input already exists. Skipping"
fi


#=============================================================
#====                       STEP 5:                       ====
#====                 Concatenating UCEs                  ====
#=============================================================

CAT_UCES=${OUTPUT}/tmp/006-concatenate-uces/
mkdir -p ${CAT_UCES}
mkdir -p ${OUTPUT}/partitioned-uces/

if [ -z "$(ls -A "${CAT_UCES}")" ]; then
	log "Concatenating UCEs..."
	UCE_LIST=$(find ${SWSC_PARSE}/*/ -type f -name "*nexus" \
		| sed 's/_right.nexus//g;s/_left.nexus//g;s/_core.nexus//g' \
		| sort | uniq)
    # create UCE dir
    for uce in $UCE_LIST; do
		uce_name=$(basename $uce)
		# check if ends with _all.nexus
		if [[ "$uce_name" == *_all.nexus ]]; then
  			warn "No flanks found for ${uce_name/_all.nexus}. It will be left as a whole"
			cp $uce ${OUTPUT}/partitioned-uces/${uce_name/_all}
		else
			mkdir -p ${CAT_UCES}/$uce_name
			cp ${uce}_left.nexus \
				${uce}_core.nexus \
				${uce}_right.nexus \
				${CAT_UCES}/$uce_name
		fi
    done
	DONEmsg
else
	warn "UCEs already prepared for concatenation. Skipping"
fi

#=============================================================
#====                       STEP 6:                       ====
#====                 Concatenating UCEs                  ====
#=============================================================


mkdir -p "${CAT_UCES}"/NEXUS/
if [ -z "$(ls -A "${CAT_UCES}/NEXUS/")" ]; then
	log "Concatenating UCEs..."
	# (2) concatenating with PHYLUCE in parallel
	UCE_DIRS=$(find ${CAT_UCES} -mindepth 1 -type d)
	N_UCES=$(find ${CAT_UCES} -mindepth 1 -type d | wc -l)
	AUX=0

	for uce in $UCE_DIRS; do
		uce_name=$(basename $uce)
		AUX=$(( AUX + 1 ))
		# concatenate with phyluce in parallel
		$CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS" \
			$CONDA_PREFIX/bin/phyluce_align_concatenate_alignments \
			--alignments "$uce" --phylip --log ${OUTPUT}/tmp/ \
			--output "${CAT_UCES}"/NEXUS/"${uce_name}" > /dev/null 2>&1
		bash "${HOME_DIR}"/progress-bar.sh $AUX "$N_UCES"
	done
	$CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
	DONEmsg

	log "Writing final nexus files..."
	# fix charsets and mv to output dir
	for uce in $(find "${CAT_UCES}"/NEXUS/ -type f); do
		uce_name=$(basename $uce)
		cat $uce \
		| sed 's/charpartition combined =.*//;s/.nexus//g' \
		| sed '1i #NEXUS' \
		> ${OUTPUT}/partitioned-uces/${uce_name}
	done
	DONEmsg
else
	warn "UCEs already concatenated. Skipping"
fi

#=============================================================
#====                       STEP 7:                       ====
#====                  Remove tmp files                   ====
#=============================================================

if test -f "${OUTPUT}/tmp/phyluce_align_concatenate_alignments.log"; then
	rm -rf ${OUTPUT}/tmp/
	DONEmsg
else
	warn "Temporary files already removed. Skipping"
fi
