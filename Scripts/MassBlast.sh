#!/bin/bash

INPUTDIR=$1
OUTPUTDIR=$2

for FILE in $INPUTDIR/*
do 
	FILEPATH=$(dirname $FILE)
	FILENAME=$(basename $FILE)
	FILENAME=${FILENAME%.*}

	if [ -d $FILE ]; then
		continue
	fi

	echo "Checking $FILENAME..."

	if [ ! -e "$OUTPUTDIR/$FILENAME.known.fasta" ]; then
		echo "$FILENAME not found, running TEBlast"
		TEBlast $FILE $OUTPUTDIR
	else
		echo "$FILENAME found, skipping"
	fi
	echo -e "\n"
done

SendText "MassBlast finished"

##This version will parse all subdirectories instead
##DIRECTORYPATH=$1

##for SUBDIRPATH in $(find $DIRECTORYPATH -maxdepth 1 -type d)
##do
	##echo "For $SUBDIRPATH"
##	SUBDIR=$(basename $SUBDIRPATH)
##	echo -e "\nChecking $SUBDIRPATH/$SUBDIR"	

##	if [ ! -e "$SUBDIRPATH/$SUBDIR.blast" ]; then
##		if [ ! -e "$SUBDIRPATH/$SUBDIR.fa" ]; then
##			echo "Missing $SUBDIR.fa, skipping $SUBDIR"
##			continue
##		fi

        	##echo "Could not locate $SUBDIRPATH/*.blast"
##		echo "Running TEBlast on $SUBDIR"
##		TEBlast "$SUBDIRPATH/$SUBDIR.fa"
	
##	else
##                echo "Found $SUBDIR.blast, skipping $SUBDIR"
##	fi
##done





SendText "MassBlast2 done running"




