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

echo "phase_svs.sh"
echo "  BINDIR: $BINDIR"
echo "  PHASEDSNPS: $PHASEDSNPS"
echo "  STRUCTURALVARIANTS: $STRUCTURALVARIANTS"
echo "  LONGREADSBAM: $LONGREADSBAM"
echo "  GENOME: $GENOME"
echo "  OUT: $OUTPREFIX"
echo
echo


if [ ! -r $OUTPREFIX.hairs ]
then
  echo "extracting pacbio-hairs from snps"
  extractHAIRS --bam $LONGREADSBAM --VCF $PHASEDSNPS --out $OUTPREFIX.hairs
fi

if [ ! -r $OUTPREFIX.spliced.vcf ]
then
  echo "Splicing in phased SVs"
  $BINDIR/splicephase.pl $PHASEDSNPS $STRUCTURALVARIANTS $OUTPREFIX.hairs $OUTPREFIX.spliced.vcf $GENOME >& $OUTPREFIX.spliced.log
fi

exit


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

