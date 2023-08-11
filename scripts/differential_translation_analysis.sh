#!/bin/bash
#### Description: Prepocesses raw sequencing data for   
####              1. Ribo-seq
####                  a) adapter trimming
####                  b) ncRNA contamination
####                  c) Read length based on phasing
####              2. RNA-seq  
####                  a) adapter trimming
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.

# Make a reference to count with Plastid
./scripts/make_plastid_reference.sh

mkdir -p quantitation/transcript_cds_rpkm
./scripts/calculate_rpkm.sh

mkdir quantitation/gene_cds_counts
./scripts/calculate_gene_cds_counts.sh

Rscript ./scripts/run_deseq2.R