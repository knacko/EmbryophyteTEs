#!/bin/bash

FILE=$1
FILEPATH=$(dirname $FILE)
FILENAME=$(basename $FILE)
FILENAME=${FILENAME%.*}
OUTPUTDIR=$2

#mkdir "$OUTPUTDIR/$FILENAME"
#cd "$OUTPUTDIR/$FILENAME"

blastn -query $FILE -db ~/libs/REdat/REdatDB -outfmt '6 length qstart qend salltitles ' -evalue 0.01 -out "$OUTPUTDIR/$FILENAME.blast" -num_alignments 1
#python ~/scripts/FormatBLASTFasta $FILEPATH/$FILENAME.blast > $OUTPUTDIR/$FILENAME.known.fasta
#rm $FILEPATH/$FILENAME.blast

##sed -r 's/.*/#/' $FILE.blast > $FILE.blast.out

##sed -r 's/^(.*?)#(.*?)#(.*?)#.*$/>\2#\3\n\1/' $FILE.blast > $FILE.blast.out

##awk 'BEGIN { OFS = "\n" } { print ">"$2"#"$3, $1 }' $FILE.blast

