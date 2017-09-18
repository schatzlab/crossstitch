#!/usr/bin/perl -w
use strict;

my $WIGGLE = 100;

my $USAGE = "comparesvcalls.pl sv.summary spliced.vcf.svphase [VERBOSE]\n";

my $TRUTH = shift @ARGV or die $USAGE;
my $SNIFFLES = shift @ARGV or die $USAGE;
my $VERBOSE = shift @ARGV;

my @truesv;

open TRUE, $TRUTH or die "Cant open $TRUTH ($!)\n";
open SNIFFLES, $SNIFFLES or die "Cant open $SNIFFLES ($!)\n";


## Load the true calls
###############################################################################

print "Loading truth ($TRUTH)\n";
while (<TRUE>)
{
  chomp;
  my ($info, $start, $chr2, $end, $type, $len) = split /\s+/, $_;
  
  my ($file, $chr1) = split /:/, $info;
  my $het = (split /\./, $file)[1];

  my $v;
  $v->{start}  = $start;
  $v->{end}    = $end;
  $v->{chr}    = $chr1;
  $v->{het}    = $het;
  $v->{len}    = $len;
  $v->{type}   = $type;
  $v->{found}  = 0;
  $v->{nearby} = 0;

  if ($VERBOSE) { print "$chr1\t$start\t$end\t$len\t$type\t$het\n"; }
  push @truesv, $v;
}

my $truecnt = scalar @truesv;
print "Loaded $truecnt true svs\n";




## Process the sniffle calls once to figure out the phase match
###############################################################################

my $header = <SNIFFLES>; # skip header

my $cnt = 0;

my %phasematch;
$phasematch{hetA}->{"0|1"} = 0;
$phasematch{hetA}->{"1|0"} = 0;
$phasematch{hetA}->{"1|1"} = 0;
$phasematch{hetB}->{"0|1"} = 0;
$phasematch{hetB}->{"1|0"} = 0;
$phasematch{hetB}->{"1|1"} = 0;

while (<SNIFFLES>)
{
  chomp;
  $cnt++;

  my ($posinfo, $type, $svlen, $seqlen, $bar1, $numreads, $hap1reads, $hap2reads, $bar2, $het, $hap1r, $bar3, $newgenotype) = split /\s+/, $_;
  my ($chr, $pos, $ogenotype) = split /:/, $posinfo;

  $type = substr($type, 1, length($type)-2);
 # print "=$cnt\t$chr\t$pos\t$svlen\t$type\t$newgenotype\n";

  foreach my $t (@truesv)
  {
    next if ($chr ne $t->{chr});
    my $distance = abs($pos - $t->{start});

    if (($distance < $WIGGLE) && ($type eq $t->{type}))
    {
      my $tchr   = $t->{chr};
      my $tstart = $t->{start};
      my $tlen   = $t->{len};
      my $thet   = $t->{het};
      my $ttype  = $t->{type};

      if ($thet ne "homAB") 
      {
        $phasematch{$thet}->{$newgenotype}++;
      }

  #    print " \t$tchr\t$tstart\t$tlen\t$ttype\t$thet\n";
    }
  }
  #print "\n";
}


## Print matching statistics
###############################################################################

my %truephase;
$truephase{"homAB"} = "1|1";

print "== printing phase matching statistics\n";
foreach my $thet (sort keys %phasematch)
{
  print "$thet: ";
  my $bestgt = "???";
  my $bestcnt = -1;

  foreach my $gt (sort keys %{$phasematch{$thet}})
  {
    my $cnt = $phasematch{$thet}->{$gt};
    print "\t$gt\t$cnt";
    if ($cnt > $bestcnt) { $bestcnt = $cnt; $bestgt = $gt; }
  }

  print "\t=>\t$bestgt $bestcnt\n";
  $truephase{$thet} = $bestgt;
}



## Now process the sniffles calls to confirm matches
###############################################################################
seek SNIFFLES, 0, 0;

print "== processing sniffles calls with phase matching\n";
$cnt=0;

$header = <SNIFFLES>;

while (<SNIFFLES>)
{
  chomp;
  $cnt++;

  my ($posinfo, $type, $svlen, $seqlen, $bar1, $numreads, $hap1reads, $hap2reads, $bar2, $het, $hap1r, $bar3, $newgenotype) = split /\s+/, $_;
  my ($chr, $pos, $ogenotype) = split /:/, $posinfo;

  $type = substr($type, 1, length($type)-2);

  my $matched = 0;
  my $summary = "=$cnt\t$chr\t$pos\t$svlen\t$type\t$newgenotype\t|\t$hap1reads\t$hap2reads\n";

  foreach my $t (@truesv)
  {
    next if ($chr ne $t->{chr});
    my $distance = abs($pos - $t->{start});

    if (($distance < $WIGGLE) && ($type eq $t->{type}))
    {
      my $tchr   = $t->{chr};
      my $tstart = $t->{start};
      my $tlen   = $t->{len};
      my $thet   = $t->{het};
      my $ttype  = $t->{type};

      my $status = "***";
      $t->{nearby}++;

      if ($truephase{$thet} eq $newgenotype)
      {
        $status = ":-)";
        $t->{found}++;
        $matched++;
      }

      $summary .= " \t$tchr\t$tstart\t$tlen\t$ttype\t$thet\t$status\n";
    }
  }

  if (($VERBOSE) || ($matched != 1))
  {
    print $summary;
  }
}




## Report any unmatched calls
###############################################################################

print "\n";
print "==Checking for any unmatched true calls\n";

foreach my $t (@truesv)
{
  if ($t->{found} == 0)
  {
      my $tchr    = $t->{chr};
      my $tstart  = $t->{start};
      my $tlen    = $t->{len};
      my $thet    = $t->{het};
      my $ttype   = $t->{type};
      my $tnearby = $t->{nearby};

      print "$tchr\t$tstart\t$tlen\t$ttype\t$thet\t$tnearby\n";
  }
}



