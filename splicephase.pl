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

my @vcfheader;
my $vcfdata;

my $vcfheaderlines = 0;
my $vcfdatalines = 0;

while (<PHASEDVCF>)
{
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

my $hairslines = 0;
while (<READHAIRS>)
{
  $hairslines++;
}

print "Loaded $hairslines hairs records\n";

my $sniffleslines = 0;
while (<SNIFFLESVCF>)
{
  $sniffleslines++;
}

print "Loaded $sniffleslines sniffles variants\n";


print OUTVCF "not done";
