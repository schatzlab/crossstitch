#!/bin/bash
set -e
min_args=1

if [ $# -eq $min_args ]; then
	 if [ -d $1 ]; then
                echo "file does not exists"a
	else
		grep -E '^#|	0\/1|	1\/0|	1\/1|	0\|1|	1\|0|	1\|1' $1 > $1'.filtered.vcf'
	 fi
else
        echo "filteres out any non homozygrous or heterozygous snps."
        echo "Required: File path of vcf file"
fi

