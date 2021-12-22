#!/bin/bash
set -e
VERSION=0.1

Usage() {
echo -e "
-----------------------------------------------------------
CURE: an automated and parallel pipeline for UCE curation
-----------------------------------------------------------
by Vinícius Franceshini-Santos & Felipe Freitas
version "$VERSION"

\e[4mUsage\e[0m:
  CURE Command ...

\e[4mCommands\e[0m:
  GeneRegion              Runs UCE curation by gene region, such as introns and exons
  UCERegion               Runs UCE curation by UCE region, such as core and flanks

"
exit
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

# "Bye" message function (green text)
BYEmsg() {
	echo -e "\033[1;32m[ CURE v$VERSION | $(date +%Y-%m-%d" "%H:%M:%S) ] ALL DONE!\033[0m"
	echo -e "
\033[1;32mThis is a goodbye message.
Probably will have citation, mail contact, and more.
Bye\033[0m"
}
# Check deps
#check_deps

tmp=$(realpath "$0")
HOME_DIR=${tmp%/*}

# Print usage if no arg
if [ $# -eq 0 ]; then Usage; fi

# Parse commands and arguments
COMMAND=$1
ARGS=$(echo "${@:2}")

case $COMMAND in
                -h | --help )
                Usage; exit
                ;;
                GeneRegion )
                ${HOME_DIR}/scripts/GeneRegion.sh --version $VERSION $args
                BYEmsg
                ;;
                UCERegion )
                ${HOME_DIR}/scripts/UCERegion.sh $args
                BYEmsg
                ;;
                * )
                error_exit "No command named $COMMAND."
                ;;
esac
