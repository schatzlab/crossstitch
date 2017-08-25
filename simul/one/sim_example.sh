#!/bin/bash

## Make sure HAPCUT, NGMLR, SURVIVOR-LRSIM, SNIFFLES, SAMTOOLS are in your path
## also conda install sbt freebayes

#set -xv

SNPDIST=100
BASE=base.fa
PARAM=simul.param
BASE=base.fa
THREADS=10

mkdir -p data


if [ ! -r data/mutA.fasta ]
then
  echo "Simulating variants"
  SURVIVOR-LRSIM 0 $BASE $PARAM 0 data/mut $SNPDIST
fi

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
  echo "reformating reads"
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
  sniffles -m data/pbAll.bam -v data/pbAll.sniffles.vcf --cluster --genotype -n $READSTOPHASE
fi


if [ ! -r data/illA.1.fq ]
then
  samtools faidx data/mutA.fasta
  numpairs=`awk '{print $2*30/200}' data/mutA.fasta.fai`
  echo "simulating $numpairs pairs for hapB"
  mason_simulator -ir data/mutA.fasta -n $numpairs -o data/illA.1.fq -or data/illA.2.fq --num-threads $THREADS
fi

if [ ! -r data/illB.1.fq ]
then
  samtools faidx data/mutB.fasta
  numpairs=`awk '{print $2*30/200}' data/mutB.fasta.fai`
  echo "simulating $numpairs pairs for hapB"
  mason_simulator -ir data/mutB.fasta -n $numpairs -o data/illB.1.fq -or data/illB.2.fq --num-threads $THREADS
fi

if [ ! -r data/illAll.1.fq ]
then
  echo "Reformat illumina reads"
 cat data/illA.1.fq data/illB.1.fq | paste - - - - | awk '{c++; print "@illread"c"/1"; print $2; print "+"; print $4}' > data/illAll.1.fq
 cat data/illA.2.fq data/illB.2.fq | paste - - - - | awk '{c++; print "@illread"c"/2"; print $2; print "+"; print $4}' > data/illAll.2.fq
fi

if [ ! -r base.fa.amb ]
then 
  echo "Indexing genome for bwa"
  bwa index base.fa
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


