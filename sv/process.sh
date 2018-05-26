y=$1
BINDIR=$2
OUTDIR=$3
bamFile=$4
fastaFile=$5
echo 'flagsflags: ' $1 ' ' $2 ' ' $3 ' ' $4 ' ' $5
numReads=`wc -l < $y`
l=${#y}
ol=${#OUTDIR}
ty=${y:ol+9:l}
c="${ty#*.}"
c="${c#*.}"
x=${y:ol+9:l-4-ol-9-${#c}-1}
echo 'Processing insertion '$c:$x

# Produce an empty file to indicate this insertion has started being processed
touch $OUTDIR/cert/$x'.txt.'$c'.cert'

# Produce a fastq file of the reads supporting this insertion
samtools view $bamFile  $c:$(($x - 10000))-$(($x + 10000)) | grep -w -f $OUTDIR/inserts/$x.txt.$c > $OUTDIR/inserts/extracted_$x.sam.$c
java -cp "${BINDIR}" ReadInsertExtraction $OUTDIR/inserts/extracted_$x.sam.$c $x > $OUTDIR/inserts/extracted_$x.fa.$c
samtools view -H $bamFile > $OUTDIR/inserts/mine_$x.bam.$c
cat $OUTDIR/inserts/extracted_$x.sam.$c >> $OUTDIR/inserts/mine_$x.bam.$c
samtools bam2fq $OUTDIR/inserts/mine_$x.bam.$c > $OUTDIR/inserts/extracted2_$x.fa.$c

# Convert from fastq to fasta format
java -cp "${BINDIR}" FastaFileFixer $OUTDIR/inserts/extracted2_$x.fa.$c $OUTDIR/inserts/fixed_$x.fa.$c
mv $OUTDIR/inserts/fixed_$x.fa.$c $OUTDIR/inserts/extracted2_$x.fa.$c

# Put the reads in a format readable by FalconSense
java -cp "${BINDIR}" FalconFormatter $OUTDIR/inserts/extracted2_$x.fa.$c $OUTDIR/inserts/extracted_falcon_$x.fa.$c
java -cp "${BINDIR}" FalconFormatter $OUTDIR/inserts/extracted_$x.fa.$c $OUTDIR/inserts/insertions_extracted_falcon_$x.fa.$c

# Run Falconsense using default parameters
quarterReads=`expr $numReads / 8`
$BINDIR/falcon_sense --min_idt 0.7 --min_len 500 --max_read_len 40896 --min_ovl_len 250 --min_cov $quarterReads --n_core 2 > $OUTDIR/falconsense_output/$x.correctedReads.fasta.$c < $OUTDIR/inserts/extracted_falcon_$x.fa.$c
$BINDIR/falcon_sense --min_idt 0.7 --min_len 500 --max_read_len 40896 --min_ovl_len 250 --min_cov 2 --n_core 2 > $OUTDIR/falconsense_output/$x.correctedReads.insertions.fasta.$c < $OUTDIR/inserts/insertions_extracted_falcon_$x.fa.$c

# Reformat the file to have the entire sequences on one line for downstream processing
java -cp "${BINDIR}" FastaFileFixer2 $OUTDIR/falconsense_output/$x.correctedReads.fasta.$c $OUTDIR/falconsense_output/$x.correctedReads.fixed.fasta.$c
mv $OUTDIR/falconsense_output/$x.correctedReads.fixed.fasta.$c $OUTDIR/falconsense_output/$x.correctedReads.fasta.$c
java -cp "${BINDIR}" FastaFileFixer2 $OUTDIR/falconsense_output/$x.correctedReads.insertions.fasta.$c $OUTDIR/falconsense_output/$x.correctedReads.insertions.fixed.fasta.$c
mv $OUTDIR/falconsense_output/$x.correctedReads.insertions.fixed.fasta.$c $OUTDIR/falconsense_output/$x.correctedReads.insertions.fasta.$c

offset=100000
left=$(($x - $offset))
if [ "$left" -lt 1 ]
then
    left=1
fi
echo $c:$left-$(($x + $offset))
echo $fastaFile
samtools faidx $fastaFile $c:$left-$(($x + $offset)) > $OUTDIR/samples/"$x".fa.$c

ngmlr -t 4 -r $OUTDIR/samples/"$x".fa.$c -q $OUTDIR/falconsense_output/"$x".correctedReads.fasta.$c -o $OUTDIR/results/"$x"_all.sam.$c
ngmlr -t 4 -r $OUTDIR/samples/"$x".fa.$c -q $OUTDIR/falconsense_output/"$x".correctedReads.insertions.fasta.$c -o $OUTDIR/results/"$x"_all_insertions.sam.$c
echo 'aligned reads'
oldinsert=`java -cp "${BINDIR}" BestInsertFinder2 $OUTDIR/results/"$x"_all.sam.$c $x $offset 'SEQ'`
insert=`java -cp "${BINDIR}" BestInsertFinder2 $OUTDIR/results/"$x"_all_insertions.sam.$c $x $offset 'SEQ'`
oldpos=`java -cp "${BINDIR}" BestInsertFinder2 $OUTDIR/results/"$x"_all.sam.$c $x $offset 'POS'`
pos=`java -cp "${BINDIR}" BestInsertFinder2 $OUTDIR/results/"$x"_all_insertions.sam.$c $x $offset 'POS'`

echo $x >> $OUTDIR/a.txt
echo 'insert: '$insert >> $OUTDIR/a.txt
echo 'oldinsert: '$oldinsert >> $OUTDIR/a.txt

echo '>insert_'$c':'$x > $OUTDIR/seqs/$x.$c.fa
echo '>insert_'$c':'$x > $OUTDIR/seqs/$x.$c.pos
if [ "$insert" != '' ]
then
    echo $insert >> $OUTDIR/seqs/$x.$c.fa
    echo $pos >> $OUTDIR/seqs/$x.$c.pos
else
    echo $oldinsert >> $OUTDIR/seqs/$x.$c.fa
    echo $oldpos >> $OUTDIR/seqs/$x.$c.pos
fi

rm $OUTDIR/samples/"$x".*

# Produce a file to indicate this insertion has finished being processed
touch $OUTDIR/cert/$x'.txt.'$c'.cert.done'
