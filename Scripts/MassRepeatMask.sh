#!/bin/bash

#Contains the fasta files for the genome(s)
INPUTDIR=$1

#The output directory that will contain folders for each fasta file processed
OUTPUTDIR=$2

#Percentage of sequences to analyze
MAXSEQ=10

#Max bases to analyze
MAXBP=10000000

for FILE in $INPUTDIR/*
do

	if [ -e $INPUTDIR/STOP ]; then
		echo "Stop flag found"
		rm $INPUTDIR/STOP
		break
	fi

        FILEPATH=$(dirname $FILE)
        FILENAME=$(basename $FILE)
        SPECIES=${FILENAME%.*}

        if [ -d $FILE ]; then
                continue
        fi

        echo "Checking $SPECIES..."

        if [ ! -d "$OUTPUTDIR/$SPECIES" ]; then
        	
		SendText "Started Masking $SPECIES"
	
		mkdir "$OUTPUTDIR/$SPECIES"
		cd "$OUTPUTDIR/$SPECIES"

		NUMSEQS=`grep '>' -o $FILE | wc -l`
		
		# Create the RepeatMasker output files (important files are .out and .tbl)
		~/tools/RepeatModeler/RepeatMasker/RepeatMasker -qq -pa 3 -lib ~/libs/REdat.fasta -a -dir "$OUTPUTDIR/$SPECIES" $FILE

		GENOMESIZE=`sed '4q;d' *.tbl | sed -r 's/^total length:\s*(.*?)\sbp.*\(.*$/\1/'`
		TESIZE=`sed '6q;d' *.tbl | sed -r 's/^bases masked:\s*(.*?)\sbp.*\(.*$/\1/'`
						
		perl ~/tools/RMparser/parseRM_simple.pl -genlen $GENOMESIZE -RMout *.fa.out -lib ~/libs/REdat.fasta			
												
#		perl ~/tools/RMparser/parseRM.pl -v -p -i "$OUTPUTDIR/$SPECIES/$FILENAME.out"
#		perl ~/tools/RMparser/parseRM.pl -v -p -l 100,1 -g $FILE -i "$OUTPUTDIR/$SPECIES/$FILENAME.out" -r ~/libs/REdat.fasta

		cp */*.all-repeats.tab ../$SPECIES.$GENOMESIZE.$TESIZE.tab
		
	else
            echo "$SPECIES directory found, skipping"
        fi
done

echo "MassRepeatMask finished"
SendText "MassRepeatMask finished"
