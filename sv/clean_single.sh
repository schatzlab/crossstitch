BINDIR=$1
OUTDIR=$2
bamFile=$3
fastaFile=$4
filename=$5
for ((i=0; i<=3; i++)); do
    if [ ! -f $filename'.done' ]; then
        nfilename="${filename##*/}"
        nfilename=${nfilename::${#nfilename}-5}
        echo 'Rerunning (parallel): '$nfilename
        timeout 2m "${BINDIR}"/process.sh $OUTDIR'/inserts/'$nfilename $BINDIR $OUTDIR $bamFile $fastaFile
    fi
done
