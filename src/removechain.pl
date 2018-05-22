#!/usr/bin/perl -w

use strict;

my $USAGE = "removechain.pl orig.chain chrA chrB chrC ... chrN > new.chain\n";

my $chainfile = shift @ARGV or die $USAGE;

my %toremove;

foreach my $chr (@ARGV)
{
  $toremove{$chr} = 1;
}

my $removecnt = scalar keys %toremove;

print STDERR "Removing $removecnt chromosomes from $chainfile: " . join(" ", keys %toremove) . "\n";

open CHAIN, "$chainfile" or die "Cant open $chainfile ($!)\n";

my $DOPRINT = 1;

while (<CHAIN>)
{
  if (/^chain/) 
  {
    my @fields = split /\s+/, $_;
    $DOPRINT = 1;
    if (exists $toremove{$fields[2]})
    { 
      $DOPRINT = 0;  
      print STDERR "Skipping chain for $fields[2]\n";
    }
  }

  print if ($DOPRINT);
}
