#!/bin/bash
#### Description: Make coverage tracks for the manuscript figures
####             1. Merged Data for Figure 2B
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie 

### 1 merged transcriptome tracks for new annotations 

# # 1. Merged Data for Figure 2B
# ## Demonstrate the coverage at transcript-level
# ## a) CHX full coverage
# ## b) CHX p-site offset
# ## c) CHX a-site offset
# ## d) HARR-ND p-site

# conda activate plastid
source ~/miniconda3/etc/profile.d/conda.sh
conda activate plastid

## sort and index transcriptome aligned BAMs
for seqtype in riboseq_harr riboseq_chx riboseq_nd
do
 if ! [ -f  data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam.bai" ]; then
    samtools sort \
    data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.bam" \
    -o data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam"

    samtools index \
    data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam"
fi
done

mkdir -p results/alignment_tracks/merged
merged_dir=results/alignment_tracks/merged

# for seqtype in riboseq_harr riboseq_nd riboseq_harr riboseq_chx
# do
#     make_wiggle \
#     --count_files data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam" \
#     --countfile_format BAM \
#     --fiveprime \
#     --offset 12 \
#     --output_format bedgraph \
#     -o $merged_dir/$seqtype.psite.transcriptome.nonorm.bedgraph 
    
# done

# for seqtype in riboseq_harr riboseq_nd riboseq_harr riboseq_chx
# do
#     make_wiggle \
#     --count_files data/$seqtype/mapped/merged/$seqtype.bam \
#     --countfile_format BAM \
#     --fiveprime \
#     --offset 12 \
#     --output_format bedgraph \
#     -o $merged_dir/$seqtype.psite.genome.nonorm.bedgraph 
    
# done

for seqtype in riboseq_chx
do
    make_wiggle \
    --count_files data/$seqtype/mapped/merged/$seqtype.bam \
    --countfile_format BAM \
    --fiveprime \
    --offset 0 \
    --output_format bedgraph \
    -o $merged_dir/$seqtype.full.genome.nonorm.bedgraph 
    
done

# for seqtype in riboseq_harr riboseq_nd riboseq_harr riboseq_chx
# do
#     make_wiggle \
#     --count_files data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam" \
#     --countfile_format BAM \
#     --fiveprime \
#     --offset 0 \
#     --output_format bedgraph \
#     --normalize \
#     -o $merged_dir/$seqtype.fullcov.transcriptome.bedgraph \
    
# done

conda deactivate


bamCoverage -b data/riboseq_chx/mapped/merged/riboseq_chx.bam  -o results/alignment_tracks/merged/riboseq_chx.fullcov.transcriptome.bedgraph --outFileFormat bedgraph -bs 1 -p 70 --normalizeUsing BPM 