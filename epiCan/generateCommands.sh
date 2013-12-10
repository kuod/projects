#!/bin/bash

#	$1: input file
#	$2: variant file
 
echo "Generating scripts on `uname -a` with args $@"

echo input=$1
echo variant=$2

export base="/cbio/grlab/nobackup/mouse_SARSFLU/"
export input="$base$1"
export log=`date +%s`

export variant=$2

echo $input 

