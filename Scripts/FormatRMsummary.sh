
function parse {

FILENAME=$(basename $1)
FILENAME=`echo "$FILENAME" | cut -d'.' -f1`

echo "Parsing $FILENAME"

#delete blank second line
sed '2d' -i $1

#remove hash character
sed '1s/^.//' -i $1

#add filename as first element
sed '1 ! s/^/'"$FILENAME"'\t/' -i $1

#add header for file anme
sed '1 s/^/species\t/' -i $1

}

for i in $1/*; do
  if [ `head -c 1 $i` == "#" ]; then
   parse $i
 fi
done
