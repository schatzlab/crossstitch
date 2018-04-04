y=$1
BINDIR=$2
OUTDIR=$3
bamFile=$4
fastaFile=$5
numReads=`wc -l < $y`
l=${#y}
ol=${#OUTDIR}
c="${y#*.}"
c="${c#*.}"
x=${y:ol+9:l-4-ol-9-${#c}-1}
echo 'Processing insertion '$x
# Produce a fastq file of the reads
echo $x
samtools view $bamFile  $c:$(($x - 10000))-$(($x + 10000)) | grep -w -f $OUTDIR/inserts/$x.txt.$c > $OUTDIR/inserts/extracted_$x.sam
java -cp "${BINDIR}" ReadInsertExtraction $OUTDIR/inserts/extracted_$x.sam $x > $OUTDIR/inserts/extracted_$x.fa
samtools view -H $bamFile > $OUTDIR/inserts/mine_$x.bam
cat $OUTDIR/inserts/extracted_$x.sam >> $OUTDIR/inserts/mine_$x.bam
samtools bam2fq $OUTDIR/inserts/mine_$x.bam > $OUTDIR/inserts/extracted2_$x.fa

# Convert from fastq to fasta format
java -cp "${BINDIR}" FastaFileFixer $OUTDIR/inserts/extracted2_$x.fa $OUTDIR/inserts/fixed_$x.fa
mv $OUTDIR/inserts/fixed_$x.fa $OUTDIR/inserts/extracted2_$x.fa
#java -cp "${BINDIR}" FastaFileFixer $OUTDIR/inserts/extracted_$x.fa $OUTDIR/inserts/insertions_fixed_$x.fa
#mv $OUTDIR/inserts/insertions_fixed_$x.fa $OUTDIR/inserts/extracted_$x.fa

# Put the reads in a format readable by FalconSense
java -cp "${BINDIR}" FalconFormatter $OUTDIR/inserts/extracted2_$x.fa $OUTDIR/inserts/extracted_falcon_$x.fa
java -cp "${BINDIR}" FalconFormatter $OUTDIR/inserts/extracted_$x.fa $OUTDIR/inserts/insertions_extracted_falcon_$x.fa

# Run Falconsense using default parameters
quarterReads=`expr $numReads / 8`
$BINDIR/falcon_sense --min_idt 0.7 --min_len 500 --max_read_len 40896 --min_ovl_len 250 --min_cov $quarterReads --n_core 2 > $OUTDIR/falconsense_output/$x.correctedReads.fasta < $OUTDIR/inserts/extracted_falcon_$x.fa
$BINDIR/falcon_sense --min_idt 0.7 --min_len 500 --max_read_len 40896 --min_ovl_len 250 --min_cov 2 --n_core 2 > $OUTDIR/falconsense_output/$x.correctedReads.insertions.fasta < $OUTDIR/inserts/insertions_extracted_falcon_$x.fa

# Reformat the file to have the entire sequences on one line for downstream processing
java -cp "${BINDIR}" FastaFileFixer2 $OUTDIR/falconsense_output/$x.correctedReads.fasta $OUTDIR/falconsense_output/$x.correctedReads.fixed.fasta
mv $OUTDIR/falconsense_output/$x.correctedReads.fixed.fasta $OUTDIR/falconsense_output/$x.correctedReads.fasta
java -cp "${BINDIR}" FastaFileFixer2 $OUTDIR/falconsense_output/$x.correctedReads.insertions.fasta $OUTDIR/falconsense_output/$x.correctedReads.insertions.fixed.fasta
mv $OUTDIR/falconsense_output/$x.correctedReads.insertions.fixed.fasta $OUTDIR/falconsense_output/$x.correctedReads.insertions.fasta

offset=100000
left=$(($x - $offset))
if [ "$left" -lt 0 ]
then
    left=0
fi
samtools faidx $fastaFile $c:$left-$(($x + $offset)) > $OUTDIR/samples/"$x".fa

touch $OUTDIR/seqs/$x.fa
echo '>insert_'$c':'$x > $OUTDIR/seqs/$x.fa

ngmlr -t 4 -r $OUTDIR/samples/"$x".fa -q $OUTDIR/falconsense_output/"$x".correctedReads.fasta -o $OUTDIR/results/"$x"_all.sam
ngmlr -t 4 -r $OUTDIR/samples/"$x".fa -q $OUTDIR/falconsense_output/"$x".correctedReads.insertions.fasta -o $OUTDIR/results/"$x"_all_insertions.sam

insert=''
oldinsert=`java -cp "${BINDIR}" BestInsertFinder2 $OUTDIR/results/"$x"_all.sam $x $offset`
insert=`java -cp "${BINDIR}" BestInsertFinder2 $OUTDIR/results/"$x"_all_insertions.sam $x $offset`

echo $x >> $OUTDIR/a.txt
echo 'insert: '$insert >> $OUTDIR/a.txt
echo 'oldinsert: '$oldinsert >> $OUTDIR/a.txt
if [ "$insert" != '' ]
then
    echo $insert >> $OUTDIR/seqs/$x.fa
else
    echo $oldinsert >> $OUTDIR/seqs/$x.fa
fi

rm $OUTDIR/samples/"$x".*
