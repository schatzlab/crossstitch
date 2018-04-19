BINDIR=`dirname $(readlink -f "$0")`
WORKINGDIR=`pwd`
mkdir $WORKINGDIR/testout
OUTDIR=$WORKINGDIR/testout
echo $WORKINGDIR
echo $BINDIR
echo $OUTDIR
rm $BINDIR/*.class
rm -r $OUTDIR/*
rm -r $OUTDIR

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

if [ -z "${vcfFile}" ] || [ -z "${bamFile}" ] || [ -z "${fastaFile}" ] || [ -z "${outputFile}" ]; then
    usage
fi

BLASR=/home-3/mkirsche\@jhu.edu/scratch/miniconda3/envs/blasr/bin/blasr
PYTHON=python

INSERT_BEFORE=1 # The number of characters before the insertion to include in the REF field of the new VCF file
INSERT_AFTER=0 # The number of characters after the insertion to include in the REF field of the new VCF file

if [ $USER == "mschatz1@jhu.edu" ]
then
BLASR=/work-zfs/mschatz1/mschatz/build/miniconda/bin/blasr
PYTHON=python2
fi


javac $BINDIR/*.java
# Generate lists of reads for all inserts
java -cp "${BINDIR}" ReadFinder $WORKINGDIR/"${vcfFile}" $OUTDIR/inserts

#for y in `ls $OUTDIR/inserts/*.txt.*`
#do
#    "${BINDIR}"/process.sh $y $BINDIR $OUTDIR $bamFile $fastaFile
#done
parallel --jobs 16 "${BINDIR}"/process.sh {} $BINDIR $OUTDIR $bamFile $fastaFile ::: $OUTDIR/inserts/*.txt.*

wait

"${BINDIR}"/clean.sh $BINDIR $OUTDIR $bamFile $fataFile
cat $OUTDIR/seqs/*.fa > $OUTDIR/all.seq
cat $OUTDIR/seqs/*.pos > $OUTDIR/all.pos
java -cp "${BINDIR}" VCFEditor $OUTDIR/all.seq $OUTDIR/all.pos $WORKINGDIR/$vcfFile $WORKINGDIR/$fastaFile $WORKINGDIR/$outputFile $INSERT_BEFORE $INSERT_AFTER

rm $WORKINGDIR/hs*.log


