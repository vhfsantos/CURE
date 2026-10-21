#!/bin/bash
set -e -o pipefail
VERSION="$2"

trap 'exit 130' SIGINT
trap 'exit 143' SIGTERM
trap 'exit 129' SIGHUP

CURRENT_STAGE=""
CURRENT_STAGE_DIRS=()
CURRENT_STAGE_LOG=""
CURRENT_STAGE_JOB_LOGS=""

collect_stage_logs() {
	[ -n "$CURRENT_STAGE_JOB_LOGS" ] || return 0
	[ -d "$CURRENT_STAGE_JOB_LOGS" ] || return 0
	local job_file
	while IFS= read -r -d '' job_file; do
		printf '\n===== %s =====\n' "$(basename "$job_file" .stdout)" >> "$CURRENT_STAGE_LOG"
		cat -- "$job_file" >> "$CURRENT_STAGE_LOG"
	done < <(find "$CURRENT_STAGE_JOB_LOGS" -maxdepth 1 -type f -name '*.stdout' -print0 | sort -zV)
	rm -rf -- "$CURRENT_STAGE_JOB_LOGS"
	CURRENT_STAGE_JOB_LOGS=""
}

run_logged_command() {
	local job_file=$1
	shift
	"$@" > "$job_file" 2>&1
}
export -f run_logged_command

# An interrupted stage is never reused: preserve completed stages and logs.
cleanup_on_exit() {
	status=$?
	trap - EXIT
	if [ "$status" -ne 0 ] && [ -n "$CURRENT_STAGE" ]; then
		# Wait for this run's queued jobs before removing their working files.
		if [ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/bin/sem" ]; then
			"$CONDA_PREFIX/bin/sem" --id "$$" --wait >/dev/null 2>&1 || true
		fi
		collect_stage_logs
		warn "Stage $CURRENT_STAGE failed; removing its temporary files. Logs are kept in $LOGDIR"
		for stage_dir in "${CURRENT_STAGE_DIRS[@]}"; do
			case "$stage_dir" in
				"$OUTPUT"/tmp/*) rm -rf -- "$stage_dir" ;;
			esac
		done
		rm -f -- "$OUTPUT/tmp/.stage-$CURRENT_STAGE.complete"
	fi
	exit "$status"
}
trap cleanup_on_exit EXIT

start_stage() {
	CURRENT_STAGE="$1"
	shift
	CURRENT_STAGE_LOG="$LOGDIR/stages/$(printf '%02d' "$CURRENT_STAGE")-$(stage_name "$CURRENT_STAGE").log"
	: > "$CURRENT_STAGE_LOG"
	CURRENT_STAGE_JOB_LOGS=""
	CURRENT_STAGE_DIRS=("$@")
	for stage_dir in "${CURRENT_STAGE_DIRS[@]}"; do
		case "$stage_dir" in
			"$OUTPUT"/tmp/*) rm -rf -- "$stage_dir"; mkdir -p -- "$stage_dir" ;;
			*) error_exit "Unsafe temporary directory: $stage_dir" ;;
		esac
	done
}

complete_stage() {
	collect_stage_logs
	touch "$OUTPUT/tmp/.stage-$CURRENT_STAGE.complete"
	CURRENT_STAGE=""
	CURRENT_STAGE_DIRS=()
	CURRENT_STAGE_LOG=""
}

stage_name() {
	case "$1" in
		0) echo prepare ;; 1) echo subgroup-concatenation ;;
		2) echo swsc ;; 3) echo parse ;;
		4) echo full-concatenation ;; 5) echo prepare-loci ;;
		6) echo partitioned-loci ;;
	esac
}

fail() {
	>&2 echo "ERROR: $1"
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
	local message="[ CURE v$VERSION | $(date +%Y-%m-%d" "%H:%M:%S) ] $1"
	echo "$message"
	[ -z "$CURRENT_STAGE_LOG" ] || echo "$message" >> "$CURRENT_STAGE_LOG"
}
# Warning message function (yellow text)
warn() {
	local message="[ CURE v$VERSION | $(date +%Y-%m-%d" "%H:%M:%S) ] $1"
	echo -e "\033[1;33m${message}\033[0m"
	[ -z "$CURRENT_STAGE_LOG" ] || echo "$message" >> "$CURRENT_STAGE_LOG"
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

[ -d "$NEXUS_DIR" ] || error_exit "NEXUS input directory does not exist: $NEXUS_DIR"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || error_exit "--threads must be a positive integer."
OUTPUT=$(realpath -m -- "$OUTPUT")
NEXUS_DIR=$(realpath -- "$NEXUS_DIR")
[ "$OUTPUT" != / ] || error_exit "The output directory cannot be /."
[ "$OUTPUT" != "$NEXUS_DIR" ] || error_exit "The output directory must differ from the input directory."
mkdir -p "$OUTPUT/tmp"
LOGDIR="$OUTPUT/logfiles/run-$(date +%Y%m%dT%H%M%S)-$$"
mkdir -p "$LOGDIR/stages"
exec > >(tee -a "$LOGDIR/pipeline.log") 2>&1


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

NEXUSCOPY="$OUTPUT/tmp/000-nexus-copy"
subgroups_dir="$OUTPUT/tmp/001-subgroups"
mapfile -d '' -t INPUT_UCES < <(find "$NEXUS_DIR" -maxdepth 1 -type f -name '*.nexus' -print0 | sort -z)
N_UCES=${#INPUT_UCES[@]}
[ "$N_UCES" -gt 0 ] || fail "No .nexus files found in $NEXUS_DIR"
n_subgroups=$THREADS
if [ "$n_subgroups" -gt "$N_UCES" ]; then n_subgroups=$N_UCES; fi

if [ ! -f "$OUTPUT/tmp/.stage-0.complete" ]; then
	start_stage 0 "$NEXUSCOPY" "$subgroups_dir"
	log "Preparing input data..."
	log "Splitting $N_UCES UCEs into $n_subgroups subgroups"
	cp -- "${INPUT_UCES[@]}" "$NEXUSCOPY/"
	for ((i=1; i<=n_subgroups; i++)); do mkdir -p "$subgroups_dir/$i"; done
	for ((i=0; i<N_UCES; i++)); do
		group=$((i % n_subgroups + 1))
		mv -- "$NEXUSCOPY/$(basename "${INPUT_UCES[$i]}")" "$subgroups_dir/$group/"
	done
	complete_stage
        DONEmsg
else
	warn "Input preparation already done. Skipping"
fi

#=============================================================
#====                       STEP 1:                       ====
#====               Concatenate UCEs for SWSC             ====
#=============================================================

SUBGROUPS_CAT="$OUTPUT/tmp/002-subgroups-cat"

if [ ! -f "$OUTPUT/tmp/.stage-1.complete" ]; then
	stage_jobs="$OUTPUT/tmp/002-subgroup-logs"
	start_stage 1 "$SUBGROUPS_CAT" "$stage_jobs"
	CURRENT_STAGE_JOB_LOGS="$stage_jobs"
	log "Concatenating UCEs subgroups for SWSC..."
        # concatenate with phyluce in parallel
        for sg in $(seq 1 $n_subgroups); do
			mkdir -p "$stage_jobs/subgroup-$sg"
                $CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS" \
				run_logged_command "$stage_jobs/subgroup-$sg.stdout" \
                        $CONDA_PREFIX/bin/phyluce_align_concatenate_alignments \
                        --alignments "${subgroups_dir}/$sg" \
	                        --nexus --log "$stage_jobs/subgroup-$sg" \
	                        --output "${SUBGROUPS_CAT}/$sg"
                bash "${HOME_DIR}"/progress-bar.sh $sg "$n_subgroups"
        done
        $CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
		collect_stage_logs
		for sg in $(seq 1 "$n_subgroups"); do
			[ -s "$SUBGROUPS_CAT/$sg/$sg.nexus" ] || fail "Missing concatenated alignment for subgroup $sg"
		done
		complete_stage
        DONEmsg
else
	warn "UCEs concatenation already done. Skipping"
fi


#=============================================================
#====                       STEP 2:                       ====
#====                      Run SWSC                       ====
#=============================================================

SWSC="$OUTPUT/tmp/003-swsc"

if [ ! -f "$OUTPUT/tmp/.stage-2.complete" ]; then
	stage_jobs="$OUTPUT/tmp/003-swsc-logs"
	start_stage 2 "$SWSC" "$stage_jobs"
	CURRENT_STAGE_JOB_LOGS="$stage_jobs"
	log "Running SWSC..."
        # run SWSC in parallel
        for sg in $(seq 1 $n_subgroups); do
                # (1) fix uce names in .nexus files
                sed -i 's/uce-/uce_/g' "${SUBGROUPS_CAT}/$sg/$sg.nexus"

                # (2) run SWSC
	                $CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS" \
				run_logged_command "$stage_jobs/subgroup-$sg.stdout" \
                        $CONDA_PREFIX/bin/python $SWSC_PATH \
                        $( realpath ${SUBGROUPS_CAT}/${sg}/${sg}.nexus ) \
	                        $( realpath $SWSC )
                bash "${HOME_DIR}"/progress-bar.sh $sg "$n_subgroups"
        done
        $CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
		collect_stage_logs
		for sg in $(seq 1 "$n_subgroups"); do
			[ -s "$SWSC/$sg.nexus_entropy_partition_finder.cfg" ] || fail "Missing SWSC output for subgroup $sg"
		done
		complete_stage
        DONEmsg
else
	warn "SWSC already run. Skipping"
fi

#=============================================================
#====                       STEP 3:                       ====
#====                   Parsing Results                   ====
#=============================================================

SWSC_PARSE="$OUTPUT/tmp/004-swsc-parse"

# function to parse results
SWSCParser(){
	subgroup=$1
	SWSC=$2
	SUBGROUPS_CAT=$3
	SWSC_PARSE=$4
	OUTPUT=$5
	PARSE_LOGDIR=$6
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
		--output ${SWSC_PARSE}/${subgroup}/ --log-path "$PARSE_LOGDIR" \
		--output-format nexus
}

# Need to export this function to call it with Parallel
export -f SWSCParser

if [ ! -f "$OUTPUT/tmp/.stage-3.complete" ]; then
	stage_jobs="$OUTPUT/tmp/004-parse-logs"
	start_stage 3 "$SWSC_PARSE" "$stage_jobs"
	CURRENT_STAGE_JOB_LOGS="$stage_jobs"
	log "Parsing results..."
        for sg in $(seq 1 $n_subgroups); do
		# calls function in parallel
			mkdir -p "$stage_jobs/subgroup-$sg"
                $CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS" \
				run_logged_command "$stage_jobs/subgroup-$sg.stdout" \
                        SWSCParser $sg ${SWSC} ${SUBGROUPS_CAT} ${SWSC_PARSE} \
							${OUTPUT} "$stage_jobs/subgroup-$sg"
                bash "${HOME_DIR}"/progress-bar.sh $sg "$n_subgroups"
        done
        $CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
		collect_stage_logs
		for sg in $(seq 1 "$n_subgroups"); do
			[ -s "$SWSC_PARSE/$sg.nexus" ] || fail "Missing parsed alignment for subgroup $sg"
			[ -n "$(find "$SWSC_PARSE/$sg" -maxdepth 1 -type f -name '*.nexus' -print -quit)" ] \
				|| fail "Missing split loci for subgroup $sg"
		done
		complete_stage
        DONEmsg
else
	warn "SWSC results already parsed. Skipping"
fi

#=============================================================
#====                       STEP 4:                       ====
#====        Generating concat files and PF2 input        ====
#=============================================================

UCES_CAT="$OUTPUT/tmp/005-concatenated-uces"


if [ ! -f "$OUTPUT/tmp/.stage-4.complete" ]; then
	start_stage 4 "$UCES_CAT"
	log "Generating concatenated files and PF2 input..."
	# (1) Place all subgroups nexus to a same dir
	for sg in $(seq 1 $n_subgroups); do
			cp "${SWSC_PARSE}/${sg}/"* "$UCES_CAT/"
	done

	# (2) Concatenate
	$CONDA_PREFIX/bin/phyluce_align_concatenate_alignments \
		--alignments "${UCES_CAT}" \
			--phylip --log "$UCES_CAT" \
			--output "${UCES_CAT}/concatenated-uces" \
			>> "$CURRENT_STAGE_LOG" 2>&1

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
	   # 3rd sed: remove 'charset ' from beggining of lines;
	   # 4rd sed: remove ".nexus" from partition names
	cat "${UCES_CAT}"/concatenated-uces/concatenated-uces.charsets \
		| sed '/charpartition combined/,$d' \
		| sed 's/begin sets;//g' \
		| sed 's/^charset //g' | tail -n +3 \
		| sed 's/.nexus//g' \
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
		[ -s "$OUTPUT/concatenated-uces/concatenated-uces-swscen.phylip" ] || fail "Missing concatenated alignment"
		[ -s "$OUTPUT/concatenated-uces/concatenated-uces-swscen.charsets" ] || fail "Missing concatenated charsets"
		complete_stage

else
	warn "PF2 input already exists. Skipping"
fi


#=============================================================
#====                       STEP 5:                       ====
#====                 Concatenating UCEs                  ====
#=============================================================

CAT_UCES="$OUTPUT/tmp/006-concatenate-uces"
mkdir -p "$OUTPUT/partitioned-uces"

if [ ! -f "$OUTPUT/tmp/.stage-5.complete" ]; then
	start_stage 5 "$CAT_UCES"
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
		complete_stage
	
else
	warn "UCEs already prepared for concatenation. Skipping"
fi

#=============================================================
#====                       STEP 6:                       ====
#====                 Concatenating UCEs                  ====
#=============================================================


if [ ! -f "$OUTPUT/tmp/.stage-6.complete" ]; then
	stage_jobs="$OUTPUT/tmp/007-partitioned-logs"
	start_stage 6 "$CAT_UCES/NEXUS" "$stage_jobs"
	CURRENT_STAGE_JOB_LOGS="$stage_jobs"
	# (2) concatenating with PHYLUCE in parallel
		UCE_DIRS=$(find "$CAT_UCES" -mindepth 1 -maxdepth 1 -type d ! -name NEXUS)
		N_UCES=$(find "$CAT_UCES" -mindepth 1 -maxdepth 1 -type d ! -name NEXUS | wc -l)
	AUX=0

	for uce in $UCE_DIRS; do
		uce_name=$(basename $uce)
		AUX=$(( AUX + 1 ))
			mkdir -p "$stage_jobs/$uce_name"
		# concatenate with phyluce in parallel
		$CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs "$THREADS" \
			run_logged_command "$stage_jobs/$uce_name.stdout" \
			$CONDA_PREFIX/bin/phyluce_align_concatenate_alignments \
				--alignments "$uce" --phylip --log "$stage_jobs/$uce_name" \
				--output "${CAT_UCES}"/NEXUS/"${uce_name}"
		bash "${HOME_DIR}"/progress-bar.sh $AUX "$N_UCES"
	done
	$CONDA_PREFIX/bin/sem --will-cite --id $$ --wait
	collect_stage_logs
		for uce in $UCE_DIRS; do
			uce_name=$(basename "$uce")
			[ -s "$CAT_UCES/NEXUS/$uce_name/$uce_name.phylip" ] || fail "Missing PHYLIP output for $uce_name"
		done
	DONEmsg
	
	log "Writing final files..."
		# mv phylip to output dir
	for uce in $(find "${CAT_UCES}"/NEXUS/ -name "*.phylip" -type f ); do
	uce_name=$(basename $uce)
		mv $uce ${OUTPUT}/partitioned-uces/${uce_name}
	done
	
	# fix charsets and mv them to output dir
	for uce in $(find "${CAT_UCES}"/NEXUS/ -name "*.charsets" -type f ); do
		uce_name=$(basename $uce)
		cat $uce \
		| sed 's/charpartition combined =.*//;s/.nexus//g' \
		| sed '1i #NEXUS' \
		> ${OUTPUT}/partitioned-uces/${uce_name}
	done
		complete_stage
	DONEmsg
	
else
	warn "UCEs already concatenated. Skipping"
fi

#=============================================================
#====                       STEP 7:                       ====
#====                  Remove tmp files                   ====
#=============================================================
	log "Removing temp files..."
rm -rf -- "$OUTPUT/tmp"
DONEmsg
echo "If you make use of the UCERegion strategy in your research, please cite the papers describing the strategy and SWSC-EN:

Freitas et al. (2021). Partitioned Gene-Tree Analyses and Gene-Based Topology Testing Help Resolve
Incongruence in a Phylogenomic Study of Host-Specialist Bees (Apidae: Eucerinae).
Molecular Biology and Evolution 38(3):1090-1100. https://doi.org/10.1093/molbev/msaa277

Tagliacollo and Lanfear (2018). Estimating Improved Partitioning Schemes for
Ultraconserved Elements. Molecular Biology and Evolution 35(7):1798-1811.
https://doi.org/10.1093/molbev/msy069"
