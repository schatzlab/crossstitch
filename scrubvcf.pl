#!/usr/bin/perl -w
use strict;

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

while (<>)
{
  if ($_ =~ /^#/) { print $_; }
  else
  {
    $all++;
    my @fields = split /\s+/, $_;
    my $chr = $fields[0];
    my $type = $fields[4];

    if (exists $chrom{$chr} && exists $types{$type})
    {
      $reported++;
      $types{$type}++;
      print $_;
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
