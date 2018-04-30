BINDIR=$1
OUTDIR=$2
bamFile=$3
fastaFile=$4
parallel --timeout 150 --jobs 4 "${BINDIR}"/clean_single.sh $BINDIR $OUTDIR $bamFile $fastaFile {} ::: $OUTDIR/cert/*.cert
