#!/bin/bash

INPUTDIR=$1

for FILE in $INPUTDIR/*
do
        FILEPATH=$(dirname $FILE)
        FILENAME=$(basename $FILE)
        FILENAME=${FILENAME%.*}
	FILEDIR=$FILENAME
	FILENAME+=".fa.out"

        echo "Checking $FILENAME..."

        if [ -f "$INPUTDIR/$FILEDIR/$FILENAME" ]; then
                echo "Running parseRM on $FILEDIR"
		##mkdir "$OUTPUTDIR/$FILENAME"
		cd "$INPUTDIR/$FILEDIR"

		perl ~/tools/RMparser/parseRM.pl -v -p -l 100,1 -g 117500000 -i $FILENAME -r ~/libs/REdat.fas
        else
                echo "$FILENAME not found, skipping"
        fi
        echo -e "\n"
done

SendText "MassRMparser finished"

