BINDIR=$1
OUTDIR=$2
bamFile=$3
fastaFile=$4
parallel --timeout 30 --jobs 16 "${BINDIR}"/clean_single.sh $BINDIR $OUTDIR $bamFile $fastaFile {} ::: $OUTDIR/inserts/*.txt.*
