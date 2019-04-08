set -e

if [ $# -ne 4 ] && [ $# -ne 5 ]
then
  echo "USAGE: refineinsertions.sh structural_variants.vcf long_reads.bam genome.fa outfile [numjobs]"
  echo ""
  echo "Details:"
  echo "  structural_variants.vcf:           VCF file of structural variants identified using Sniffles"
  echo "  long_reads.bam:                    BAM file of long reads aligned with NGMLR"
  echo "  genome.fa:                         Reference genome used"
  echo "  outfile:                           Name of output VCF file"
  echo "  numjobs:                           The number of jobs to run when refining insertions"
  exit
fi

jobs=16

if [ $# -eq 5 ]
then
    jobs=$5
fi

BINDIR=`dirname $(readlink -f "$0")`
$BINDIR/../sv/go.sh -v $1 -b $2 -f $3 -o $4 -j $jobs
