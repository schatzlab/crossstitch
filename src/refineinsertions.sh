set -e

if [ $# -ne 4 ]
then
  echo "USAGE: refineinsertions.sh structural_variants.vcf long_reads.bam genome.fa outfile"
  echo ""
  echo "Details:"
  echo "  structural_variants.vcf:           VCF file of structural variants identified using Sniffles"
  echo "  long_reads.bam:                    BAM file of long reads aligned with NGMLR"
  echo "  genome.fa:                         Reference genome used"
  echo "  outfile:                           Name of output VCF file"
  exit
fi

BINDIR=`dirname $(readlink -f "$0")`
$BINDIR/../sv/go.sh -v $1 -b $2 -f $3 -o $4
