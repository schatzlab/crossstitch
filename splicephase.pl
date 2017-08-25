#!/usr/bin/perl -w
use strict;

my $USAGE = "splicephase.pl phased.vcf sniffles.vcf loadreads.hairs spliced.vcf\n";

my $PHASEDVCFFILE   = shift or die $USAGE;
my $SNIFFLESVCFFILE = shift or die $USAGE;
my $READHAIRSFILE   = shift or die $USAGE;
my $OUTVCFFILE      = shift or die $USAGE;


open PHASEDVCF,   $PHASEDVCFFILE   or die "Cant open $PHASEDVCFFILE ($!)\n";
open SNIFFLESVCF, $SNIFFLESVCFFILE or die "Cant open $SNIFFLESVCFFILE ($!)\n";
open READHAIRS,   $READHAIRSFILE   or die "Cant open $READHAIRSFILE ($!)\n";
open OUTVCF,       "> $OUTVCFFILE" or die "Cant open $OUTVCFFILE ($!)\n";





## Load the phased VCF file
###############################################################################

my @vcfheader;
my $vcfdata;

my $vcfheaderlines = 0;
my $vcfdatalines = 0;

while (<PHASEDVCF>)
{
  chomp;
  if (/^#/) 
  {
    push @vcfheader, $_;
    $vcfheaderlines++;
  }
  else
  {
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = split /\s+/, $_; 
    my $v;
    $v->{chrom}  = $chrom;
    $v->{pos}    = $pos;
    $v->{id}     = $id;
    $v->{ref}    = $ref;
    $v->{alt}    = $alt;
    $v->{qual}   = $qual;
    $v->{filter} = $filter;
    $v->{info}   = $info;
    $v->{format} = $format;
    $v->{sample} = $sample;

    $vcfdata->{$chrom}->{$pos} = $v;
    $vcfdatalines++;
  }
}

print "Loaded $vcfheaderlines header lines and $vcfdatalines variants\n";


## Load Sniffles SV calls
###############################################################################

my $readstophase;
my $snifflesvariants;
my $sniffleslines = 0;
while (<SNIFFLESVCF>)
{
  if (/^#/) 
  {
    ## nothing to do
  }
  else
  {
    chomp;
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = split /\s+/, $_; 
    my $v;
    $v->{chrom}  = $chrom;
    $v->{pos}    = $pos;
    $v->{id}     = $id;
    $v->{ref}    = $ref;
    $v->{alt}    = $alt;
    $v->{qual}   = $qual;
    $v->{filter} = $filter;
    $v->{info}   = $info;
    $v->{format} = $format;
    $v->{sample} = $sample;
    $v->{reads}  = [];

    my @infofields = split /;/, $info;
    foreach my $f (@infofields)
    {
      if ($f =~ /^RNAMES/)
      {
         $f = substr($f, 7);
         my @reads = split /,/, $f;
         foreach my $r (@reads)
         {
           $readstophase->{$r}->{num}++;
           $readstophase->{$r}->{sv}->{$chrom}->{$pos}->{phase} = -1;
           push @{$v->{reads}}, $r;
         }
      }
    }

    $snifflesvariants->{$chrom}->{$pos} = $v;
    $sniffleslines++;
  }
}

my $readcount = scalar keys %$readstophase;
print "Loaded $sniffleslines sniffles variants involving $readcount reads\n";

foreach my $r (sort keys %$readstophase)
{
	my $num = $readstophase->{$r}->{num};
	print "  $r $num\n";
}


## Determine the read phase information
###############################################################################

my $hairslines = 0;
while (<READHAIRS>)
{
  $hairslines++;
}

print "Loaded $hairslines hairs records\n";


print OUTVCF "not done";
