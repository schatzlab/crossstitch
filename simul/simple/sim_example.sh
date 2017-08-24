#!/bin/sh

## Make sure HAPCUT, NGMLR, SURVIVOR-LRSIM, SNIFFLES, SAMTOOLS are in your path

SNP_DIST=100
BASE=base.fa
SURVIVOR-LRSIM 0 base.fa simul.param 0 mut $SNPDIST
