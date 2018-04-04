# CrossStitch: Hybrid Phasing and Personal Genome Construction

CrossStitch creates personalized reference-quality diploid genomes without de novo assembly. The basic idea is rather than trying to assemble a genome from scratch, it will leverage a reference genome as a baseline, and then update it with any SNPs, indels, or structural variations present in your sample. For the best results, the data requirements are similar to a de novo assembly: Illumina-based data for SNPs and Indels, Long Read data for structural variants, and Phasing data such as 10X Linked Reads and/or HiC data. However the CrossStitch process is much less demanding, produces more accurate results, and the process is much more predictable. The output will be a phased VCF file with all variants (SNPs, Indels, and SVs) as well as a phased personalized diploid genome including 2 copies of each chromosome with the variants inserted at the correct locations.

[Slides describing CrossStitch](http://schatz-lab.org/presentations/2018/2018.PAG.DiploidGenomes.ShortLongLinkedReads.pdf)



## Installation

```

## CrossStitch requires extractHairs from HapCut2
$ git clone https://github.com/vibansal/HapCUT2
$ cd HapCUT2
$ make

$ git clone --recursive https://github.com/schatzlab/crossstitch.git
$ cd vcf2diploid
$ make

## After this you will need to add HapCUT2/build/extractHAIRS to your path or edit src/crossstitch.sh with the correct path

```


## Running CrossStitch

Currently only human genomes are supported for diploid genome construction.

```
$ crossstich.sh phased_snps.vcf unphased_structural_variants.vcf long_reads.bam genome.fa outputprefix gender refine
 
Details:
  phased_snps.vcf:                   VCF file of phased SNP and indel variants. Recommend LongRanger (10X only) or HapCUT2 (HiC and/or 10X)
  unphased_structural_variants.vcf:  VCF file of structural variants identified using Sniffles
  long_reads.bam:                    BAM file of long reads aligned with NGMLR
  genome.fa:                         Reference genome used
  outputprefix:                      Prefix for output files
  gender:                            "male" or "female", used to ensure sex chromosomes are correctly used
  refine:                            optionally refine structural variant calls with local assembly (1=refine, 0=skip)
```


## Output files

### Main Files: 
```
*.hap1.fa.gz: Haplotype 1 chromosome fasta sequences 
*.hap2.fa.gz: Haplotype 2 chromosome fasts sequences 
*.spliced.scrubbed.vcf.gz: Finalized set of phased small and structural variants 
```

### Annotation files: 
```
*.map:   liftover file to relate coordinates on the personalized assembly to the reference (such as GRCh38) 
*.chain: liftover file to relate coordinates on the personalized assembly to the reference (such as GRCh38) 
```



