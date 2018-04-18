BINDIR=$1
OUTDIR=$2
bamFile=$3
fastaFile=$4

for filename in $OUTDIR/cert/*.cert; do
    for ((i=0; i<=3; i++)); do
        if [ ! -f $filename.done ]; then
            nfilename="${filename##*/}"
            nfilename=${nfilename::-4}
            "${BINDIR}"/process.sh {} $BINDIR $OUTDIR $bamFile $fastaFile $nfilename
        fi
    done
done
