#!/bin/bash

#set -xv
set -e

if [ $# -ne 5 ]
then
  echo "USAGE: phase_svs.sh phased_snps.vcf unphased_structural_variants.vcf long_reads.bam genome.fa outputprefix"
  exit
fi

BINDIR=`dirname $(readlink -f "$0")`
PHASEDSNPS=$1
STRUCTURALVARIANTS=$2
LONGREADSBAM=$3
GENOME=$4
OUTPREFIX=$5
VCF2DIPLOIDJAR=/work-zfs/mschatz1/mschatz/build/vcf2diploid/vcf2diploid.jar

echo "phase_svs.sh"
echo "  BINDIR: $BINDIR"
echo "  PHASEDSNPS: $PHASEDSNPS"
echo "  STRUCTURALVARIANTS: $STRUCTURALVARIANTS"
echo "  LONGREADSBAM: $LONGREADSBAM"
echo "  GENOME: $GENOME"
echo "  OUT: $OUTPREFIX"
echo
echo

GENOME=`readlink -f $GENOME`


if [ ! -r $OUTPREFIX.hairs ]
then
  echo "extracting pacbio-hairs from phased snps"
  extractHAIRS --bam $LONGREADSBAM --VCF $PHASEDSNPS --out $OUTPREFIX.hairs
fi

if [ ! -r $OUTPREFIX.scrubbed.vcf ]
then
  echo "Scrubbing SV calls"
  ($BINDIR/scrubvcf.pl $STRUCTURALVARIANTS > $OUTPREFIX.scrubbed.vcf) >& $OUTPREFIX.scrubbed.log
fi

if [ ! -r $OUTPREFIX.spliced.vcf ]
then
  echo "Splicing in phased SVs"
  $BINDIR/splicephase.pl $PHASEDSNPS $OUTPREFIX.scrubbed.vcf $OUTPREFIX.hairs $OUTPREFIX.spliced.vcf $GENOME >& $OUTPREFIX.spliced.log
fi

if [ ! -r $OUTPREFIX.alleleseq ]
then
  mkdir -p $OUTPREFIX.alleleseq
  cd $OUTPREFIX.alleleseq

  ln -s ../$OUTPREFIX.spliced.vcf

  VCFID=`head -5000 $OUTPREFIX.spliced.vcf | grep '#CHROM' | awk '{print $10}'`

  echo "constructing diploid sequence with SNPs and SVs"
  java -Xmx400000m -jar $VCF2DIPLOIDJAR -id $VCFID -pass -chr $GENOME -vcf $OUTPREFIX.spliced.vcf >& vcf2diploid.log
  cd ..
fi

if [ ! -r $OUTPREFIX.alleleseq.raw.tgz ]
then
  echo "tarring up alleleseq"
  tar czvf $OUTPREFIX.alleleseq.raw.tgz $OUTPREFIX.alleleseq
fi

