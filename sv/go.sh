BINDIR=`dirname $(readlink -f "$0")`
WORKINGDIR=`pwd`
mkdir $WORKINGDIR/csout
OUTDIR=$WORKINGDIR/csout
echo $WORKINGDIR
echo $BINDIR
echo $OUTDIR
rm -r $OUTDIR/*
rm -r $OUTDIR
export _JAVA_OPTIONS="-XX:ParallelGCThreads=16"
javac $BINDIR/*.java
mkdir $OUTDIR
mkdir $OUTDIR/alignments
mkdir $OUTDIR/falconsense_output
mkdir $OUTDIR/assembly
mkdir $OUTDIR/inserts
mkdir $OUTDIR/results
mkdir $OUTDIR/seqs
mkdir $OUTDIR/samples
mkdir $OUTDIR/cert

usage() { echo "Usage: $0 -v <vcfFile> -b <bamFile> -f <fastaFile> -o <outputFile>" 1>&2; exit 1; }

while getopts v:b:f:o: option
do
    case "${option}"
    in
    v) vcfFile=${OPTARG};;
    b) bamFile=${OPTARG};;
    f) fastaFile=${OPTARG};;
    o) outputFile=$OPTARG;;
 esac
done
echo $fastaFile
if [ -z "${vcfFile}" ] || [ -z "${bamFile}" ] || [ -z "${fastaFile}" ] || [ -z "${outputFile}" ]; then
    usage
fi

vcfPath=$WORKINGDIR/$vcfFile
if [[ $vcfFile == /* ]]; then
    vcfPath=$vcfFile
fi

fastaPath=$WORKINGDIR/$fastaFile
if [[ $fastaFile == /* ]]; then
    fastaPath=$fastaFile
fi

INSERT_BEFORE=1 # The number of characters before the insertion to include in the REF field of the new VCF file
INSERT_AFTER=0 # The number of characters after the insertion to include in the REF field of the new VCF file

# Generate lists of reads for all insertions
java -cp "${BINDIR}" ReadFinder $vcfPath $OUTDIR/inserts

# Process all insertions in parallel
parallel --timeout 500 --jobs 16 "${BINDIR}"/process.sh {} $BINDIR $OUTDIR $bamFile $fastaPath ::: $OUTDIR/inserts/*.txt.*

wait

"${BINDIR}"/clean_parallel.sh $BINDIR $OUTDIR $bamFile $fastaFile
cat $OUTDIR/seqs/*.fa > $OUTDIR/all.seq
cat $OUTDIR/seqs/*.pos > $OUTDIR/all.pos
java -cp "${BINDIR}" VCFEditor $OUTDIR/all.seq $OUTDIR/all.pos $vcfPath $fastaPath $WORKINGDIR/$outputFile $INSERT_BEFORE $INSERT_AFTER


