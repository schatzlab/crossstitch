#!/bin/bash

#set -xv
set -e

if [ $# -ne 7 ]
then
  echo "USAGE: crossstitch.sh phased_snps.vcf unphased_structural_variants.vcf long_reads.bam genome.fa outputprefix karyotype refine"
  echo ""
  echo "Details:"
  echo "  phased_snps.vcf:                   VCF file of phased SNP and indel variants. Recommend LongRanger (10X only) or HapCUT2 (HiC and/or 10X)"
  echo "  unphased_structural_variants.vcf:  VCF file of structural variants identified using Sniffles"
  echo "  long_reads.bam:                    BAM file of long reads aligned with NGMLR"
  echo "  genome.fa:                         Reference genome used"
  echo "  outputprefix:                      Prefix for output files"
  echo "  karyotype:                         "xy" or "xx", used to ensure sex chromosomes are correctly used"
  echo "  refine:                            optionally refine structural variant calls with local assembly (1=refine, 0=skip)"
  exit
fi

BINDIR=`cd "$(dirname "$0")" ; pwd -P`
PHASEDSNPS=$1
STRUCTURALVARIANTS=$2
LONGREADSBAM=$3
GENOME=$4
OUTPREFIX=$5
KARYOTYPE=$6
REFINE=$7

VCF2DIPLOIDJAR=$BINDIR/../vcf2diploid/vcf2diploid.jar
EXTRACTHAIRS=/work-zfs/mschatz1/mkirsche/github/HapCUT2/build/extractHAIRS
GZIP=pigz

if [ $USER = "mschatz" ]
then
  EXTRACTHAIRS=extractHAIRS
fi

echo "crossstich.sh"
echo "  BINDIR: $BINDIR"
echo "  PHASEDSNPS: $PHASEDSNPS"
echo "  STRUCTURALVARIANTS: $STRUCTURALVARIANTS"
echo "  LONGREADSBAM: $LONGREADSBAM"
echo "  GENOME: $GENOME"
echo "  KARYOTYPE: $KARYOTYPE"
echo "  REFINE: $REFINE"
echo 
echo "  OUT: $OUTPREFIX"
echo
echo

## Sanity check parameters

if [[ $KARYOTYPE != "xy" && $KARYOTYPE != "xx" ]]
then
  echo "Unknown karyotype: $KARYOTYPE (must be xy or xx)"
  exit
fi

if [ ! -r $GENOME ]
then
  echo "Cannot read genome file: $GENOME"
  exit
fi

if [ ! -r $LONGREADSBAM ]
then
  echo "Cannot read long reads bam file: $LONGREADSBAM"
  exit
fi

if [ ! -r $STRUCTURALVARIANTS ]
then
  echo "Cannot read sv vcf file: $STRUCTURALVARIANTS"
  exit
fi

if [ ! -r $PHASEDSNPS ]
then
  echo "Cannot read phased snp file: $PHASEDSNPS"
  exit
fi

## Sanity checks passed, begin analysis

GENOME="$(cd "$(dirname "$GENOME")" ; pwd -P)/$(basename $GENOME)"
AS=$OUTPREFIX.alleleseq
VCFID=`head -5000 $PHASEDSNPS | grep '#CHROM' | awk '{print $10}'`

javac $BINDIR/*.java

if [ ! -r $OUTPREFIX.hairs ]
then
  if [ ! -r $PHASEDSNPS.prehairs ]
  then
    echo "preprocessing phased snps to remove haplotype genotype calls"
    java -cp "${BINDIR}" RemoveStrayHairs $PHASEDSNPS $PHASEDSNPS.prehairs
  fi
  echo "extracting pacbio-hairs from phased snps (mbq 4)"
  $EXTRACTHAIRS --mbq 4 --bam $LONGREADSBAM --VCF $PHASEDSNPS.prehairs --out $OUTPREFIX.hairs
fi

if [[ $REFINE == "1" ]]
then 
  if [ ! -r $OUTPREFIX.refined.vcf ]
  then
    echo "Refining SVs"
    #$BINDIR/../RefineInsertions/rebuild_external.sh
    $BINDIR/../Iris/build.sh
    java -cp $BINDIR/../Iris/src Iris genome_in=$GENOME vcf_in=$STRUCTURALVARIANTS reads_in=$LONGREADSBAM vcf_out=$OUTPREFIX.refined.vcf threads=12
  fi
else
  if [ ! -r $OUTPREFIX.refined.vcf ]
  then
    echo "Skip SV refinement"
    cp $STRUCTURALVARIANTS $OUTPREFIX.refined.vcf
  fi
fi

if [ ! -r $OUTPREFIX.scrubbed.vcf ]
then
  echo "Scrubbing SV calls"
  (java -cp $BINDIR RemoveInvalidVariants $OUTPREFIX.refined.vcf $OUTPREFIX.scrubbed.vcf) >& $OUTPREFIX.scrubbed.log
fi

if [ ! -r $OUTPREFIX.spliced.vcf ]
then
  echo "Splicing in phased SVs"
  java -cp $BINDIR PhaseSVs $PHASEDSNPS $OUTPREFIX.scrubbed.vcf $OUTPREFIX.hairs $OUTPREFIX.spliced.vcf $GENOME >& $OUTPREFIX.spliced.log
fi

if [ ! -r $OUTPREFIX.spliced.scrubbed.vcf ]
then
  echo "Final scrub to remove overlapping spliced variants"
  (java -cp $BINDIR RemoveInvalidVariants -o $KARYOTYPE $OUTPREFIX.spliced.vcf $OUTPREFIX.spliced.scrubbed.vcf) >& $OUTPREFIX.spliced.scrubbed.log
fi

if [ ! -r $OUTPREFIX.spliced.scrubbed.vcf.gz ]
then
  echo "compressing spliced scrubbed vcf"
  $GZIP -c $OUTPREFIX.spliced.scrubbed.vcf > $OUTPREFIX.spliced.scrubbed.vcf.gz
fi

if [ ! -r $AS ]
then
  mkdir -p $AS
  pushd $AS

  echo "constructing diploid sequence with SNPs and SVs"
  java -Xmx400000m -jar $VCF2DIPLOIDJAR -id $VCFID -pass -chr $GENOME -vcf ../$OUTPREFIX.spliced.scrubbed.vcf >& vcf2diploid.log
  popd
fi

if [ ! -r $AS.raw.tgz ]
then
  echo "tarring up alleleseq"
  tar czvf $AS.raw.tgz $AS
fi

if [ ! -r $AS/raw/ ]
then
  echo "renaming files"

  pushd $AS

  mkdir -p raw
  mkdir -p raw/attic

  for i in `ls *_$VCFID*`
  do 
    echo " $i"
    j=`echo $i | sed "s/_$VCFID//"`
    mv $i $OUTPREFIX.$j
  done

  for i in `/bin/ls *_paternal*`
  do 
    j=`echo $i | sed s/_paternal/.hap1/`
    mv $i $j
  done

  for i in `/bin/ls *_maternal*`
  do 
    j=`echo $i | sed s/_maternal/.hap2/`
    mv $i $j
  done

  mv paternal.chain hap1.chain
  mv maternal.chain hap2.chain

  mv hap1.chain raw/$OUTPREFIX.hap1.chain
  mv hap2.chain raw/$OUTPREFIX.hap2.chain

  if [ -r $OUTPREFIX.chrM.hap2.fa ]
  then
    mv $OUTPREFIX.chrM.hap2.fa raw/attic
  fi

  if [[ $KARYOTYPE == "xy" ]]
  then
    echo "xy sample, making X and Y haploid"

    if [ -r $OUTPREFIX.chrX.hap2.fa ]
    then
      mv *chrX.hap2.fa *chrY.hap2.fa raw/attic
    fi

    cat raw/$OUTPREFIX.hap1.chain                                   | sed 's/paternal/hap1/' > $OUTPREFIX.hap1.chain
    $BINDIR/removechain.pl raw/$OUTPREFIX.hap2.chain chrM chrX chrY | sed 's/maternal/hap2/' > $OUTPREFIX.hap2.chain
  else
    echo "xx sample, stashing Y chromosome"

    if [ -r $OUTPREFIX.chrY.hap1.fa ]
    then
      mv *chrY* raw/attic
    fi

    $BINDIR/removechain.pl raw/$OUTPREFIX.hap1.chain chrY           | sed 's/paternal/hap1/' > $OUTPREFIX.hap1.chain
    $BINDIR/removechain.pl raw/$OUTPREFIX.hap2.chain chrM chrY      | sed 's/maternal/hap2/' > $OUTPREFIX.hap2.chain
  fi

  echo "fixing map files"
  for i in `/bin/ls *.map`
  do
    echo " $i"
    mv $i raw/
    sed 's/PAT/HAP1/' raw/$i | sed 's/MAT/HAP2/' > $i
  done

  echo "fixing chromosome files"
  for i in `/bin/ls *hap1.fa`
  do
    echo " $i"
    mv $i raw/
    sed 's/paternal/hap1/' raw/$i > $i
  done

  for i in `/bin/ls *hap2.fa`
  do
    echo " $i"
    mv $i raw/
    sed 's/maternal/hap2/' raw/$i > $i
  done

  popd
fi

if [ ! -r $AS/$OUTPREFIX.hap1.fa.gz ]
then
  echo "Making final diploid genome"
  cat `ls $AS/*.chr*.hap1.fa | sort -V` > $AS/$OUTPREFIX.hap1.fa &
  cat `ls $AS/*.chr*.hap2.fa | sort -V` > $AS/$OUTPREFIX.hap2.fa &

  wait

  echo "compressing genome files"
  $GZIP -c $AS/$OUTPREFIX.hap1.fa > $AS/$OUTPREFIX.hap1.fa.gz
  $GZIP -c $AS/$OUTPREFIX.hap2.fa > $AS/$OUTPREFIX.hap2.fa.gz
fi
