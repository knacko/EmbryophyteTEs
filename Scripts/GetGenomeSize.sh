#!/bin/bash
for FILE in "$@"
do
	FILENAME=$(basename $FILE)
	FILENAME=${FILENAME%.*}
	echo -e "$FILENAME: \c"
	grep -v ">" $1 | wc | awk '{print $3-$1}' 
done
