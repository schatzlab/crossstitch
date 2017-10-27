#!/usr/bin/perl -w
use strict;

my $STRIP_OVERLAP = 0;

if (scalar @ARGV > 0)
{
  if ($ARGV[0] eq "-h")
  {
    print "USAGE: scrubvcf.pl [-o] orig.vcf > new.vcf\n";
    print "\n";
    print " By default, scrub variants that arent on chr1-22, chrX, chrY or are not INS, DEL, INV\n";
    print " In -o mode, remove any overlapping variants or variants not on chr1-22, chrX, chrY\n";
    print "             but DONT scrub any variant types\n";

    exit(0);
  }
  elsif ($ARGV[0] eq "-o")
  {
    print STDERR "stripping overlapping calls\n";
    $STRIP_OVERLAP = 1;
    shift @ARGV;
  }
}

## Only report variants on chr1 - chr22, chrX, chrY
my %chrom;
for (my $i = 1; $i <= 22; $i++)
{
  $chrom{"chr$i"} = 0;
}

$chrom{"chrX"} = 0;
$chrom{"chrY"} = 0;

## Types to report
my %types;
$types{"<INS>"} = 0;
$types{"<DEL>"} = 0;
$types{"<INV>"} = 0;

my $all = 0;
my $reported = 0;

my $lastchr = undef;

my $lastpos0 = -1; my $lastref0 = ""; my $lastalt0 = "";
my $lastpos1 = -1; my $lastref1 = ""; my $lastalt1 = "";

while (<>)
{
  if ($_ =~ /^#/) { print $_; }
  else
  {
    $all++;

    my @fields   = split /\s+/, $_;
    my $chr      = $fields[0];
    my $pos      = $fields[1];
    my $id       = $fields[2];
    my $ref      = $fields[3];
    my $alt      = $fields[4];
    my $qual     = $fields[5];
    my $filter   = $fields[6];
    my $genotype = (split (/:/, $fields[9]))[0];

    if ($STRIP_OVERLAP)
    {
      if ($filter ne "PASS")
      {
        $types{"FAIL_PASS"}++;
        next;
      }

      if (!exists $chrom{$chr})
      {
        $types{"FAIL_CHROM"}++;
        next;
      }

      my $printvar = 1;

      my $type = "SUB";
      if    (length($ref) < length($alt)) { $type = "INS"; }
      elsif (length($ref) > length($alt)) { $type = "DEL"; }

      my $hap = 2; ## all unphased genotypes and 1|1 are on both haplotypes
      if    ($genotype eq "0|1") { $hap = 1; }
      elsif ($genotype eq "1|0") { $hap = 0; }

      if ((defined $lastchr) && ($chr eq $lastchr))
      {
        if (($hap == 0) || ($hap == 2)) { if ($pos < $lastpos0) 
        { $printvar = 0; $types{"FAIL_$hap"}++; print STDERR " overlap detected on hap $hap (lastpos: $lastpos0 $lastref0 $lastalt0) at $_"; } }

        if (($hap == 1) || ($hap == 2)) { if ($pos < $lastpos1) 
        { $printvar = 0; $types{"FAIL_$hap"}++; print STDERR " overlap detected on hap $hap (lastpos: $lastpos1 $lastref1 $lastalt1) at $_"; } }
      } 

      if ($printvar)
      {
        $reported++;
        $types{$type}++;
        print $_;

        if ((!defined $lastchr) || ($lastchr ne $chr))
        {
          # reset the chromosome
          $lastchr = $chr;
          $lastpos0 = -1; $lastref0 = ""; $lastalt0 = "";
          $lastpos1 = -1; $lastref1 = ""; $lastalt1 = "";
        }

        my $newpos = $pos + length($ref) - 1;

        if (($hap == 0) || ($hap == 2)) { $lastpos0 = $newpos; $lastref0 = $ref; $lastalt0 = $alt; }
        if (($hap == 1) || ($hap == 2)) { $lastpos1 = $newpos; $lastref1 = $ref; $lastalt1 = $alt; }
      }
    }
    else
    {
      if (exists $chrom{$chr} && exists $types{$alt})
      {
        $reported++;
        $types{$alt}++;
        print $_;
      }
    }
  }
}

print STDERR "Reported $reported of $all variants:";
foreach my $t (sort keys %types)
{
  my $n = $types{$t};
  print STDERR " $t $n";
}
print STDERR "\n";
