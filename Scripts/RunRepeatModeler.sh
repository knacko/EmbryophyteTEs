#!/bin/bash

FILE=$1
OUTPUTDIR=$2

FILENAME=$(basename $FILE)
FILENAME=${FILENAME%.*}

mkdir "$OUTPUTDIR/$FILENAME"
cd "$OUTPUTDIR/$FILENAME"

BuildDatabase -name $FILENAME -engine ncbi $FILE
time RepeatModeler -engine ncbi -pa 2 -database $FILENAME >& run.out

SendText "RepeatModeler finished"
