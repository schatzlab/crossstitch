#!/bin/bash

set -xv

## Make sure HAPCUT, NGMLR, SURVIVOR-LRSIM, SNIFFLES, SAMTOOLS are in your path

SNPDIST=100
BASE=base.fa
PARAM=simul.param
BASE=base.fa

mkdir -p data

SURVIVOR-LRSIM 0 $BASE $PARAM 0 data/mut $SNPDIST
