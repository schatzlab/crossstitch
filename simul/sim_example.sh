#!/bin/bash

## Make sure HAPCUT, NGMLR, SURVIVOR-LRSIM, SNIFFLES, SAMTOOLS are in your path
## also conda install sbt freebayes

#set -xv

SNPDIST=1000
BASE=base.fa
PARAM=simul.param
BASE=base.fa
THREADS=10

if [ ! -r $PARAM ]
then
  echo "Run me from one of the test directories!"
  exit
fi

mkdir -p data

if [ ! -r base.fa.amb ]
then 
  echo "Indexing genome for bwa"
  bwa index base.fa
fi

if [ ! -r data/mutA.fasta ]
then
  echo "Simulating variants"
  SURVIVOR-LRSIM 0 $BASE $PARAM 0 data/mut $SNPDIST
fi

if [ ! -r data/svs.summary ]
then
  grep -v SNP data/mut*.bed  | sort -nk2 | awk '{print $0, $4-$2+1}' > data/svs.summary
fi

## Simulate PacBio Reads
###############################################################################

if [ ! -r data/pbA.fa ]
then
  echo "Simulating PacBio reads for hapA"
  SURVIVOR 2 simul data/mutA.fasta ~/build/SURVIVOR/HG002_Pac_error_profile_bwa.txt 30 data/pbA.fa
fi

if [ ! -r data/pbB.fa ]
then
  echo "Simulating PacBio reads for hapB"
  SURVIVOR 2 simul data/mutB.fasta ~/build/SURVIVOR/HG002_Pac_error_profile_bwa.txt 30 data/pbB.fa
fi

if [ ! -r data/pbAll.fa ]
then
  echo "reformating pacbio reads"
  awk '{if (/^>/){c++; print ">pbread"c;}else{print}}' data/pbA.fa data/pbB.fa > data/pbAll.fa
fi

if [ ! -r data/pbAll.sam ]
then
  echo "Aligning PacBio reads and sorting"
  ngmlr -t $THREADS -r base.fa -q data/pbAll.fa -o data/pbAll.sam
  samtools view -b data/pbAll.sam -o data/pbAll.unsorted.bam
  samtools sort data/pbAll.unsorted.bam -o data/pbAll.bam
  samtools index data/pbAll.bam
  rm -f data/pbAll.unsorted.bam
fi

if [ ! -r data/pbAll.sniffles.vcf ]
then
  echo "Calling variants with Sniffles"
  READSTOPHASE=1000
  /work-zfs/mschatz1/mschatz/build/Sniffles/bin/sniffles-core-1.0.6/sniffles -m data/pbAll.bam -v data/pbAll.sniffles.vcf --cluster --genotype --report_seq -n $READSTOPHASE
fi


## Simulate Illumina Reads
###############################################################################


if [ ! -r data/illA.1.fq ]
then
  samtools faidx data/mutA.fasta
  numpairs=`awk '{print $2*50/200}' data/mutA.fasta.fai`
  echo "simulating $numpairs pairs for hapA"
  mason_simulator -ir data/mutA.fasta -n $numpairs -o data/illA.1.fq -or data/illA.2.fq --num-threads $THREADS
fi

if [ ! -r data/illB.1.fq ]
then
  samtools faidx data/mutB.fasta
  numpairs=`awk '{print $2*50/200}' data/mutB.fasta.fai`
  echo "simulating $numpairs pairs for hapB"
  mason_simulator -ir data/mutB.fasta -n $numpairs -o data/illB.1.fq -or data/illB.2.fq --num-threads $THREADS
fi

if [ ! -r data/illAll.1.fq ]
then
  echo "Reformat illumina reads"
 cat data/illA.1.fq data/illB.1.fq | paste - - - - | awk '{c++; print "@illread"c"/1"; print $2; print "+"; print $4}' > data/illAll.1.fq
 cat data/illA.2.fq data/illB.2.fq | paste - - - - | awk '{c++; print "@illread"c"/2"; print $2; print "+"; print $4}' > data/illAll.2.fq
fi


if [ ! -r data/illAll.bam ]
then
  echo "Align illumina reads and sort"
  bwa mem base.fa data/illAll.1.fq data/illAll.2.fq -t $THREADS > data/illAll.sam
  samtools view -b data/illAll.sam -o data/illAll.unsorted.bam
  samtools sort data/illAll.unsorted.bam -o data/illAll.bam
  samtools index data/illAll.bam
  rm -f data/illAll.unsorted.bam
fi

if [ ! -r data/illAll.vcf ]
then 
  echo "running freebayes"
  MIN_READS_FOR_SNP=5
  freebayes -f base.fa data/illAll.bam -C $MIN_READS_FOR_SNP > data/illAll.vcf
fi


## Phase PacBio reads using Illumina SNPs
###############################################################################

if [ ! -r data/pbAll.hairs ]
then
  echo "extracting pacbio-hairs from snps"
  extractHAIRS --bam data/pbAll.bam --VCF data/illAll.vcf --out data/pbAll.hairs
fi

if [ ! -r data/pbAll.hapcut ]
then
  echo "phase pacbio reads"
  HAPCUT2 --fragments data/pbAll.hairs --vcf data/illAll.vcf --output data/pbAll.hapcut
fi


if [ ! -r data/pbAll.phased.vcf ]
then
  echo "Making a new phased vcf file from pacbio reads"
  java -jar ~/build/fgbio/target/scala-2.12/fgbio-0.2.1-SNAPSHOT.jar HapCutToVcf -i data/pbAll.hapcut -v data/illAll.vcf -o data/pbAll.phased.vcf
fi

if [ ! -r data/pbsnps/maternal.chain ]
then
  mkdir -p data/pbsnps
  cd data/pbsnps
  ln -s ../pbAll.phased.vcf
  ln -s ../../base.fa

  echo "constructing diploid sequence with SNPs"
  java -jar /work-zfs/mschatz1/mschatz/build/vcf2diploid/vcf2diploid.jar -id unknown -chr base.fa -vcf pbAll.phased.vcf
  cd ../..
fi


if [ ! -r data/spliced.vcf ]
then
  echo "Splicing in phased SVs"
  ../../splicephase.pl data/pbAll.phased.vcf data/pbAll.sniffles.vcf data/pbAll.hairs data/spliced.vcf
fi


if [ ! -r data/spliceddiploid/maternal.chain ]
then
  mkdir -p data/spliceddiploid
  cd data/spliceddiploid
  ln -s ../spliced.vcf
  ln -s ../../base.fa

  echo "constructing diploid sequence with SNPs and SVs"
  java -jar /work-zfs/mschatz1/mschatz/build/vcf2diploid/vcf2diploid.jar -id unknown -chr base.fa -vcf spliced.vcf
  cd ../..
fi

if [ ! -r data/spliceddiploid/Am.delta ]
then
  echo "aligning diploid to truth"
  nucmer -maxmatch -D 10 data/mutA.fasta data/spliceddiploid/chr1_unknown_maternal.fa -p data/spliceddiploid/Am >& /dev/null
  nucmer -maxmatch -D 10 data/mutB.fasta data/spliceddiploid/chr1_unknown_maternal.fa -p data/spliceddiploid/Bm >& /dev/null
  nucmer -maxmatch -D 10 data/mutA.fasta data/spliceddiploid/chr1_unknown_paternal.fa -p data/spliceddiploid/Ap >& /dev/null
  nucmer -maxmatch -D 10 data/mutB.fasta data/spliceddiploid/chr1_unknown_paternal.fa -p data/spliceddiploid/Bp >& /dev/null
fi

if [ ! -r data/spliceddiploid/Am.1delta ]
then
  echo "filtering alignments"
  delta-filter -1 data/spliceddiploid/Am.delta > data/spliceddiploid/Am.1delta
  delta-filter -1 data/spliceddiploid/Bm.delta > data/spliceddiploid/Bm.1delta
  delta-filter -1 data/spliceddiploid/Ap.delta > data/spliceddiploid/Ap.1delta
  delta-filter -1 data/spliceddiploid/Bp.delta > data/spliceddiploid/Bp.1delta
fi



SHOWCOORDS=0

if [  $# -gt 0 ]
then 
  SHOWCOORDS=$1
fi

if [ $SHOWCOORDS == "1" ]
then
  echo "Am"; show-coords -rcl data/spliceddiploid/Am.1delta
  echo "Bm"; show-coords -rcl data/spliceddiploid/Bm.1delta
  echo "Ap"; show-coords -rcl data/spliceddiploid/Ap.1delta
  echo "Bp"; show-coords -rcl data/spliceddiploid/Bp.1delta
fi

