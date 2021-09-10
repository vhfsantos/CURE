#!/bin/bash

# thanks stackoverflow!
# https://stackoverflow.com/questions/238073/how-to-add-a-progress-bar-to-a-shell-script

let _progress=(${1}*100/${2}*100)/100
let _done=(${_progress})/10
let _left=10-$_done

_fill=$(printf "%${_done}s")
_empty=$(printf "%${_left}s")

printf "|${_fill// /▇}${_empty// /-}| ${_progress}%% ($1 of $2)\r"
