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
-----------------------------------------------------------
CURE: an automated and parallel pipeline for UCE curation
-----------------------------------------------------------
by Vinícius Franceshini-Santos & Felipe Freitas
version "$VERSION"

\e[4mUsage\e[0m:
 CURE UCERegion  --phyluce-nexus <nexus_dir> --output <output_dir>

\e[4mRequired arguments\e[0m:
  -p, --phyluce-nexus     Path to the directory with nexus files created with
                          phyluce
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

# Setting home dir for utils calling
tmp=$(realpath "$0")
HOME_DIR=${tmp%/*}
THREADS=2
ONLY_BY_GENE="False"
ONLY_BY_REGION="False"
FILTER="_"

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

tmp=$(realpath "$0")
HOME_DIR=${tmp%/*}


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
n_subgroups=$((`find ${NEXUS_DIR} -maxdepth 1 -type f | wc -l`/16))
echo $n_subgroups
if [ -z "$(ls -A "${NEXUSCOPY}")" ]; then
	log "Preparing input data..."
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

if [ -z "$(ls -A "${SUBGROUPS_CAT}")" ]; then
	log "Concatenating UCEs subgroups for SWSC..."
        # concatenate with phyluce in parallel
        for sg in $(seq 1 $n_subgroups); do
                $CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS" \
                        $CONDA_PREFIX/bin/phyluce_align_concatenate_alignments \
                        --alignments "${subgroups_dir}/$sg" \
                        --nexus --log ${OUTPUT}/tmp/ \
                        --output "${SUBGROUPS_CAT}/$sg" > /dev/null 2>&1
                "${HOME_DIR}"/progress-bar.sh $sg "$n_subgroups"
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

SWSC_PATH="../PFinderUCE-SWSC-EN-master/py_script/SWSCEN.py"

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
                        $( realpath $SWSC ) > /dev/null 2>&1
                "${HOME_DIR}"/progress-bar.sh $sg "$n_subgroups"
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

SWSC_PARSE=${OUTPUT}/tmp/004-swsc-parse/
mkdir -p ${SWSC_PARSE}

# assume all UCE names start with this prefix
UCE_PREFIX="uce_"

# function to parse results
SWSCParser(){
	subgroup=$1
	# (1) extract flanks coordenates for this subgroup
	   # 1st sed: adds 'charset' at the beggining of each line
	   # 2nd sed: adds 'begin sets;' as first line
	   # 3rd sed: adds 'end;' as last line
	grep $UCE_PREFIX ${SWSC}/${subgroup}.nexus_entropy_partition_finder.cfg \
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
	phyluce_align_split_concat_nexus_to_loci \
		--nexus ${SWSC_PARSE}/${subgroup}.nexus \
		--output ${SWSC_PARSE}/${subgroup}/ --log-path ${OUTPUT}/tmp/\
		--output-format nexus
}

if [ -z "$(ls -A "${SWSC_PARSE}")" ]; then
	log "Parsing results..."
        for sg in $(seq 1 $n_subgroups); do
		# calls function in parallel
                $CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS" \
                        SWSCParser $sg > /dev/null 2>&1
                "${HOME_DIR}"/progress-bar.sh $sg "$n_subgroups"
        done
        $CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
        DONEmsg
else
	warn "SWSC results already parsed. Skipping"
fi
