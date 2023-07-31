#!/bin/bash
#### Description: Prepocesses raw sequencing data for   
####              1. Ribo-seq
####                  a) adapter trimming
####                  b) ncRNA contamination
####                  c) Read length based on phasing
####              2. RNA-seq  
####                  a) adapter trimming
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# num_cores=$(getconf_NPROCESSORS_ONLN)
set -x

star_path=../bin/STAR-2.7.8a/bin/Linux_x86_64

# build the index
if ! [ -d reference_genome/star_index_riboseq ]; then
mkdir reference_genome/star_index_riboseq
$star_path/STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --sjdbOverhang 31 \
     --genomeChrBinNbits 16 \
     --genomeDir reference_genome/star_index_riboseq \
     --genomeFastaFiles reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
     --sjdbGTFfile reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf
fi

# Ribo-seq
# build indexes for contamination filtering
contam_dir=reference_genome/contaminant_sequences

# a) rRNA
if ! [ -f $contam_dir/rRNA.1.ebwt ]; then
    bowtie-build $contam_dir/rna_central_v22_rRNA.fasta $contam_dir/rRNA
fi

# b) tRNA
if ! [ -f $contam_dir/tRNA.1.ebwt ]; then
    bowtie-build $contam_dir/rna_central_v22_tRNA.fasta $contam_dir/tRNA
fi

# c) snoRNA
if ! [ -f $contam_dir/snoRNA.1.ebwt ]; then
    bowtie-build $contam_dir/rna_central_v22_snoRNA.fasta $contam_dir/snoRNA
fi

# trim and filter contamination from each sample 
# for each type of Ribo-seq data
for seqtype in riboseq_chx riboseq_harr riboseq_nd
do
   mkdir -p data/$seqtype/preproccessed_data/trimmed
   mkdir  data/$seqtype/preproccessed_data/rRNA_filter
   mkdir  data/$seqtype/preproccessed_data/tRNA_filter
   mkdir  data/$seqtype/preproccessed_data/snoRNA_filter
   mkdir data/$seqtype/preproccessed_data/complete
   mkdir -p data/$seqtype/mapped/individual

while read -ra a ;
  do
    # run cutadapt; note chx data has adapter removed by sequencing provider / -m flag used to remove 
    # short reads post trimming for HARR/ND

    if ! [ -f data/$seqtype/preproccessed_data/trimmed/${a[0]} ]; then   
      if [ $seqtype == 'riboseq_chx' ]
      then
        cutadapt  --report=full -a AGATCGGAAGAGCACACGTCT -j 50 -m 20 \
        -o data/$seqtype/preproccessed_data/trimmed/${a[0]} data/$seqtype/raw_data/${a[0]}
      else
        cutadapt  --discard-untrimmed -m 20 --report=full -a AGATCGGAAGAGCACACGTCT -j 50 --minimum-length 1 \
        -o data/$seqtype/preproccessed_data/trimmed/${a[0]} data/$seqtype/raw_data/${a[0]}
      fi
    fi
  
# Perform action here
  
    # a) filter rRNA
    if ! [ -f data/$seqtype/preproccessed_data/rRNA_filter/${a[1]}_unaligned.fq ]; then
        gunzip -c data/$seqtype/preproccessed_data/trimmed/${a[0]} | \
        bowtie -v 2 -p 50 -l 20 -norc \
        $contam_dir/rRNA \
        -q - \
        data/$seqtype/preproccessed_data/rRNA_filter/${a[1]}_aligned.fq \
        --un data/$seqtype/preproccessed_data/rRNA_filter/${a[1]}_unaligned.fq
    fi

    # b) filter snoRNA
    if ! [ -f data/$seqtype/preproccessed_data/snoRNA_filter/${a[1]}_unaligned.fq ]; then
        bowtie -v 2 -p 50 -l 20 -norc \
        $contam_dir/snoRNA \
        -q data/$seqtype/preproccessed_data/rRNA_filter/${a[1]}_unaligned.fq \
        data/$seqtype/preproccessed_data/snoRNA_filter/${a[1]}_aligned.fq \
        --un data/$seqtype/preproccessed_data/snoRNA_filter/${a[1]}_unaligned.fq
    fi

    # c) filter tRNA
    if ! [ -f data/$seqtype/preproccessed_data/tRNA_filter/${a[1]}_unaligned.fq ]; then
       bowtie -v 2 -p 50 -l 20 -norc \
        $contam_dir/tRNA \
        -q data/$seqtype/preproccessed_data/snoRNA_filter/${a[1]}_unaligned.fq \
        data/$seqtype/preproccessed_data/tRNA_filter/${a[1]}_aligned.fq \
        --un data/$seqtype/preproccessed_data/tRNA_filter/${a[1]}_unaligned.fq
    fi

    # d) remove unecessary files
    # rm data/$seqtype/preproccessed_data/*_filter/${a[1]}.*.Aligned.out.sam
    # rm data/$seqtype/preproccessed_data/*_filter/${a[1]}.*.Log.progress.out
    # rm data/$seqtype/preproccessed_data/*_filter/${a[1]}.*.SJ.out.tab
  

    # e) filter the read lengths based on phasing qc (28nt to 31nt)
    if ! [ -f data/$seqtype/preproccessed_data/complete/${a[1]}.fastq ]; then
      awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 28 && length(seq) <= 31) {print header, seq, qheader, qseq}}' \
      < data/$seqtype/preproccessed_data/tRNA_filter/${a[1]}_unaligned.fq > data/$seqtype/preproccessed_data/complete/${a[1]}.fastq
    fi
    # map the retained RPFs to the PICR reference genome

    echo "mapping to genome"

    $star_path/STAR \
    --outFilterType BySJout \
    --runThreadN 16 \
    --seedSearchStartLmaxOverLread .5 \
    --outFilterMismatchNmax 2 \
    --genomeDir reference_genome/star_index_riboseq \
    --readFilesIn data/$seqtype/preproccessed_data/complete/${a[1]}.fastq \
    --outFileNamePrefix data/$seqtype/mapped/individual/${a[1]} \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outSAMattributes All

    # rename
    mv data/$seqtype/mapped/individual/${a[1]}"Aligned.sortedByCoord.out.bam" \
    data/$seqtype/mapped/individual/${a[1]}".bam"

    # index bam
    samtools index data/$seqtype/mapped/individual/${a[1]}".bam"

  done < data/"$seqtype".txt
done

# merge and map the files for each treatment for ORF-RATER
for seqtype in riboseq_harr riboseq_nd riboseq_chx
  do
    mkdir data/$seqtype/mapped/merged

    # merge the RPF for each Ribo-seq type 
    cat data/$seqtype/preproccessed_data/complete/*.fastq > \
    data/$seqtype/mapped/merged/$seqtype.fastq

    $star_path/STAR \
    --outFilterType BySJout \
    --runThreadN 16 \
    --outFilterMismatchNmax 2 \
    --seedSearchStartLmaxOverLread .5 \
    --genomeDir reference_genome/star_index_riboseq \
    --readFilesIn data/$seqtype/mapped/merged/$seqtype.fastq \
    --outFileNamePrefix data/$seqtype/mapped/merged/$seqtype \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outSAMattributes All

    # rename
    mv data/$seqtype/mapped/merged/$seqtype"Aligned.sortedByCoord.out.bam" \
    data/$seqtype/mapped/merged/$seqtype".bam"

    # index bam
    samtools index data/$seqtype/mapped/merged/$seqtype.bam
done

# 2. RNASeq data 


mkdir -p data/rnaseq_se/preprocessed_data/complete
mkdir -p data/rnaseq_se/mapped/individual

# individual samples
while read -ra a ;
 do
     # run cutadapt
     cutadapt  -q 30 -m 20 --report=full \
     -a AGATCGGAAGAGCACACGTCT -j 32 \
     --minimum-length 1 \
     -o data/rnaseq_se/preprocessed_data/complete/${a[0]} data/rnaseq_se/raw_data/${a[0]}

    $star_path/STAR \
       --outSAMtype BAM SortedByCoordinate \
       --runThreadN 16 \
       --outFilterMismatchNmax 2 \
       --seedSearchStartLmaxOverLread .5 \
       --genomeDir reference_genome/star_index_riboseq \
       --readFilesIn  data/rnaseq_se/preprocessed_data/complete/${a[0]} \
       --outFileNamePrefix data/rnaseq_se/mapped/individual/${a[1]} \
       --outFilterMultimapNmax 1 \
       --outFilterMatchNmin 16 \
       --alignEndsType EndToEnd \
       --readFilesCommand zcat \
       --outMultimapperOrder Random \
       --outSAMattributes All

     mv data/rnaseq_se/mapped/individual/${a[1]}Aligned.sortedByCoord.out.bam data/rnaseq_se/mapped/individual/${a[1]}.bam
     samtools index data/rnaseq_se/mapped/individual/${a[1]}.bam
done < data/rnaseq_se.txt

# map the merged RNA-seq reads
mkdir -p data/rnaseq_se/mapped/merged 

zcat data/rnaseq_se/preprocessed_data/complete/* > \
data/rnaseq_se/mapped/merged/rnaseq_se.fastq 

$star_path/STAR \
    --outFilterType BySJout \
    --runThreadN 16 \
    --seedSearchStartLmaxOverLread .5 \
    --outFilterMismatchNmax 2 \
    --genomeDir reference_genome/star_index_riboseq \
    --readFilesIn data/rnaseq_se/mapped/merged/rnaseq_se.fastq \
    --outFileNamePrefix data/rnaseq_se/mapped/merged/rnaseq_se \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outSAMattributes All

    mv data/rnaseq_se/mapped/merged/rnaseq_seAligned.sortedByCoord.out.bam \
    data/rnaseq_se/mapped/merged/rnaseq_se.bam

    # index bam
    samtools index data/rnaseq_se/mapped/merged/rnaseq_se.bam
# end

# count the reads and determine read length distribution 
./scripts/fastq_read_count.sh