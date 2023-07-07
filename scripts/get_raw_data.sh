#!/bin/bash
#### Description: Downloads the RIBO-seq and RNA-seq data 
####              from ENA using FFQ
#### 
#### Written by: Clarke Lab. NIBRT - colin.clarke@nibrt.ie 

# 1. Total RNASeq (Single-end)
rna_dir=data/rnaseq_se/raw_data
mkdir -p $rna_dir

rna_samples=("SRR16796972" "SRR16796971" "SRR16796960" "SRR16796949" 
             "SRR16796946" "SRR16796945" "SRR16796944" "SRR16796943")

for sample in ${rna_samples[@]}; do
    (cd $rna_dir && ffq --ftp $sample | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs curl -O -s)
done

# 2. CHX riboseq
chx_dir=data/riboseq_chx/raw_data
mkdir -p $chx_dir

chx_samples=("SRR16796955" "SRR16796954" "SRR16796953" "SRR16796952" 
             "SRR16796951" "SRR16796950" "SRR16796948" "SRR16796947")

for sample in ${chx_samples[@]}; do
    (cd $chx_dir && ffq --ftp $sample | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs curl -O -s)
done

# 3. Harr riboseq
harr_dir=data/riboseq_harr/raw_data
mkdir -p $harr_dir

harr_samples=("SRR16796942" "SRR16796941" "SRR16796970" "SRR16796969" 
              "SRR16796968" "SRR16796967" "SRR16796966" "SRR16796965")

for sample in ${harr_samples[@]}; do
    (cd $harr_dir && ffq --ftp $sample | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs curl -O -s)
done

# 4. No drug riboseq
nd_dir=data/riboseq_nd/raw_data
mkdir -p $nd_dir

nd_samples=("SRR16796964" "SRR16796963" "SRR16796962" "SRR16796961" 
            "SRR16796959" "SRR16796958" "SRR16796957" "SRR16796956")

for sample in ${nd_samples[@]}; do
    (cd $nd_dir && ffq --ftp $sample | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs curl -O -s)
done