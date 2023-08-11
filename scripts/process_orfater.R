#!/usr/bin/env Rscript --vanilla
#### Description: Filter ORF-RATER output
####              
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# load libraries
package_list <- c("tidyverse","ORFik", "zoo")

invisible(lapply(package_list, require, character.only = TRUE, quietly=TRUE))
source("scripts/orf_filtering_functions.R")

# Annotations
## NCBI txdb
### Make NCBI reference txdbs

# make txDB
if (!file.exists("reference_genome/PICRH_ORFik.gtf.db")) {
    # make a new GTF file compatible with ORFik
    system("sed 's/XM_/XM./g; s/XR_/XR./g; s/NR_/NR./g; s/NM_/NM./g' reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf | grep -v unassigned_transcript  > reference_genome/PICRH_ORFik.gtf", intern = TRUE)

    gtf <- "reference_genome/PICRH_ORFik.gtf"
    genome <- "reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna"
    organism <- "Cricetulus griseus"

    makeTxdbFromGenome(gtf, 
        genome, 
        organism, 
        optimize = TRUE, pseudo_5UTRS_if_needed = 100)

}




### Import NCBI Feature table
suppressWarnings(
  suppressMessages(
    CriGri_PICRH_1_0_annotation <- read_delim(paste0("reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_feature_table.txt"), 
          "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
      filter(str_detect(product_accession, "XM_|NM_|NR_|XR_")) %>% 
      mutate(tid= product_accession)
    )
  )

invisible(txdb <- loadTxdb("reference_genome/PICRH_ORFik.gtf.db",
chrStyle = NULL))

# BAM files
## Genome aligned
### Ribo-seq and offset
shifts <- data.frame(
    cbind(
        fraction=c(28,29,30,31),
        offsets_start=c(-12,-12,-12,-12)
        )
        )

print(paste0("Loading genome and transcriptome BAM files for Ribo-seq and RNA-seq"))

bam_file_chx_genome <- "data/riboseq_chx/mapped/merged/riboseq_chx.bam"
genome_footprints_chx <- readBam(bam_file_chx_genome)
shiftedFootprints_chx <- shiftFootprints(genome_footprints_chx,shifts)

bam_file_harr_genome <- "data/riboseq_harr/mapped/merged/riboseq_harr.bam"
genome_footprints_harr <- readBam(bam_file_harr_genome)
shiftedFootprints_harr <- shiftFootprints(genome_footprints_harr,shifts)

bam_file_nd_genome <- "data/riboseq_nd/mapped/merged/riboseq_nd.bam"
genome_footprints_nd <- readBam(bam_file_nd_genome)
shiftedFootprints_nd <- shiftFootprints(genome_footprints_nd,shifts)

### RNA-seq
bam_file_rna <- "data/rnaseq_se/mapped/merged/rnaseq_se.bam"
genome_rna <- readBam(bam_file_rna)


## Transcriptome aligned
### Ribo-seq

bam_file_chx_transcriptome <- "data/riboseq_chx/mapped/merged/riboseq_chxAligned.toTranscriptome.out.bam"
footprints_chx_transcriptome <- readBam(bam_file_chx_transcriptome)

bam_file_harr <- "data/riboseq_harr/mapped/merged/riboseq_harrAligned.toTranscriptome.out.bam"
footprints_harr_transcriptome  <- readBam(bam_file_harr)

bam_file_nd <- "data/riboseq_nd/mapped/merged/riboseq_ndAligned.toTranscriptome.out.bam"
footprints_nd_transcriptome  <- readBam(bam_file_nd)


# ORF-RATER Predictions
## Load

orftable <- read.csv(
    file = "orfrater_analysis/rate_regression.ncbi.csv", header = T
    ) %>% 
    filter(
        orfrating >=0.5 & AAlen >= 5
        )

print(paste0(dim(orftable)[1], " ORFs with an orfrating >= 0.5 and longer than 5aa"))

# Filtering
## 1. Remove low confident ORF types
orftable <- orftable %>%
  filter(!orftype %in% c("Ciso", "Giso", "NCiso","Niso", "new_iso","internal","LOOF","truncation"))

print(paste0(dim(orftable)[1], " remain following removal of truncated, internal, LOOF and low confidence isoforms"))

## Annotate
# add annotation, select and rename columns
orftable_annotated <- orftable %>%
  left_join(CriGri_PICRH_1_0_annotation, by="tid") %>%
  dplyr::select(
    `ORF-RATER name` = orfname,
    `ORF type` = orftype,
    `ORF-RATER score` = orfrating,
    `Associated Gene symbol` = symbol, 
    `Associated Gene Name` = name, 
    `Transcript family` = tfam,
    `Transcript ID` = tid,
    `Start codon` = codon,
    `Length (AAs)` = AAlen,
    `Start Annotated?` = annot_start,
    `Stop Annotated` = annot_stop,
    `Transcript start position` = tcoord,
    `Transcript stop position` =  tstop,
    `Chromosome` = chrom,
    `Strand` = strand.x,
    `Genomic start position` = gcoord,
    `Genomic stop position` =  gstop)

# Rename ORF types
orftable_annotated <- orftable_annotated %>% 
  mutate(`ORF type` = case_when(
    `ORF type` ==  "Xiso" ~ "Isoform",
    `ORF type` ==  "Siso" ~ "Isoform",
    `ORF type` == "annotated" ~ "Annotated",
    `ORF type` == "new" ~ "New",
    `ORF type` == "upstream" ~ "Upstream",
    `ORF type` == "downstream" ~ "Downstream",
    `ORF type` == "start_overlap" ~ "Start overlap",
    `ORF type` == "stop_overlap" ~ "Stop overlap",
    `ORF type` == "extension" ~ "Extension",
  ))

## 2. TIS enrichement filtering
### Transcript RPF counts
print("Filtering ORFs by Harr enrichment at TIS")
print("Counting Harringtonine and No drug RPFs aligned to transcripts")

# Reads per transcript in transcriptome aligned BAM
harr_transcript_counts <- count_transcript_alignments(footprints_harr_transcriptome)
nd_transcript_counts <- count_transcript_alignments(footprints_nd_transcriptome)

transcript_counts <- harr_transcript_counts %>%
  dplyr::rename(nharr=counts) %>%
  left_join(nd_transcript_counts, by = "transcript") %>%
  dplyr::rename(nnd=counts) %>%
  dplyr::rename("Transcript ID"=transcript)

### TIS counts
print("Counting Harringtonine and No drug RPFs aligned to TIS")
orf_start_codons <- GRanges(
  seqnames = orftable_annotated$`Transcript ID`,
  ranges = IRanges(start=orftable_annotated$`Transcript start position`,
                   width = 4,
                   strand = "*"),
  names=orftable_annotated$`ORF-RATER name`
  )

orf_tis_counts_harr <- countOverlaps(orf_start_codons, footprints_harr_transcriptome)
orf_tis_counts_nd <- countOverlaps(orf_start_codons, footprints_nd_transcriptome)


### Rharr-Rnd calculation
orftable_annotated <- bind_cols(orftable_annotated,
                            data.frame(xharr = orf_tis_counts_harr,
                                       xnd = orf_tis_counts_nd)) %>%
  left_join(transcript_counts, by = "Transcript ID") %>%
  mutate(Rharr = (xharr/nharr)*10,
         Rnd = (xnd/nnd)*10,
         Rharr_min_Rnd=Rharr-Rnd) %>%
  mutate(tis_enriched = case_when(
  `ORF type` != "Annotated" & Rharr_min_Rnd >= 0.01 ~ "TRUE",
  `ORF type` == "Annotated" & Rharr_min_Rnd >= 0.01 ~ "TRUE",
  TRUE ~ "FALSE"
  )) %>%
  filter(tis_enriched == "TRUE") %>%
  dplyr::select(-c(tis_enriched, Rharr, Rnd))


paste0(dim(orftable_annotated)[1], " remain following TIS enirchment filtering")

## 3. FLOSS score
### Create ORF CDS GRanges
print("Calculating FLOSS score for ORFs")



# add a new annotation for each ORF
orftable_annotated <- orftable_annotated %>% 
  arrange(`Transcript ID`) %>%
  mutate(orfik_id = ave(
    `Transcript ID`, `Transcript ID`, 
    FUN = function(x) ifelse(duplicated(x), paste0(x, "_", seq_along(x)), paste0(x, "_1")))) %>%
  mutate(orfik_id =  sub("_",".", orfik_id))

### Import ORF-RATER txdb
txdb_orf <- loadTxdb("orfrater_analysis/orfrater_predictions.reference.gtf.db", chrStyle = NULL)

# Granges with ORF-Rater CDS
cds_novel_orfs <- loadRegion(txdb_orf, part = "cds")

# index to reorder
order_indices <- match(orftable_annotated$`ORF-RATER name`, names(cds_novel_orfs))

sum(orftable_annotated$`ORF-RATER name` == names(cds_novel_orfs[order_indices]))/length(orftable_annotated$`ORF-RATER name`)

cds_novel_orfs <- cds_novel_orfs[order_indices]

names(cds_novel_orfs) <- orftable_annotated$orfik_id

### Subset CDS ranges
#### Reference ORFs
reference_orfs <- orftable_annotated %>%
  filter(`ORF type` == "Annotated")

reference_cds_orfs <- cds_novel_orfs[names(cds_novel_orfs) %in% reference_orfs$orfik_id]

#### Upstream ORFs
upstream_novel_orfs <- orftable_annotated %>%
  filter(`ORF type` == "Upstream" | `ORF type` == "Start overlap") 

upstream_novel_orfs_cds <- cds_novel_orfs[names(cds_novel_orfs) %in% upstream_novel_orfs$orfik_id]

#### Other novel ORFs
other_novel_orfs <- orftable_annotated %>%
  filter(`ORF type` != "Annotated" & `ORF type` != "Upstream" & `ORF type` != "Start overlap")

other_novel_orfs_cds <- cds_novel_orfs[names(cds_novel_orfs) %in% other_novel_orfs$orfik_id]

### Calculate FLOSS using ORFik
#### Reference ORFs

orfik_chx_reference_orfs <- computeFeatures(grl = reference_cds_orfs,
                      RFP = shiftedFootprints_chx,
                      RNA = genome_rna, 
                      Gtf = txdb, 
                      faFile = "reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna", 
                      sequenceFeatures = T,
                      riboStart = 28,
                      riboStop = 31,
                      uorfFeatures = F) 

reference_orfs <- bind_cols(reference_orfs, orfik_chx_reference_orfs)

#### Upstream ORFs
orfik_chx_upstream_orfs <- computeFeatures(grl = upstream_novel_orfs_cds,
                      RFP = shiftedFootprints_chx, 
                      RNA = genome_rna, 
                      Gtf = txdb, 
                      faFile = "reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna", 
                      sequenceFeatures = T,
                      riboStart = 28,
                      riboStop = 31,
                      uorfFeatures = T) 

upstream_novel_orfs <- bind_cols(upstream_novel_orfs, orfik_chx_upstream_orfs) 

#### Other novel ORFs
orfik_chx_other_orfs <- computeFeatures(grl = other_novel_orfs_cds,
                        RFP = shiftedFootprints_chx, 
                        RNA = genome_rna, 
                        Gtf = txdb, 
                        faFile = "reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna", 
                        sequenceFeatures = T,
                        riboStart = 28,
                        riboStop = 31,
                        uorfFeatures = F)

other_novel_orfs <- bind_cols(other_novel_orfs, orfik_chx_other_orfs)

#### Merge
orftable_annotated_backup <- orftable_annotated
orftable_annotated <- orftable_annotated_backup
orftable_annotated <- bind_rows(reference_orfs, 
                                upstream_novel_orfs, 
                                other_novel_orfs)
```

## Floss Classification

### Score and count RPFs on NCBI annotated CDS
# import the annotated CDS from NCBI
cds <- loadRegion(txdb, part = "cds")

ncbi_annotated_orfs <- computeFeatures(grl = cds,
                      RFP = shiftedFootprints_chx,
                      RNA = genome_rna, 
                      Gtf = txdb, 
                      faFile = "reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna", 
                      sequenceFeatures = T,
                      riboStart = 28,
                      riboStop = 31,
                      uorfFeatures = F)

### Create model
### Classify annotated ORFs
classifier_ref <- flossCutoffs(ncbi_annotated_orfs$countRFP[ncbi_annotated_orfs$countRFP > 10],
                               ncbi_annotated_orfs$floss[ncbi_annotated_orfs$countRFP > 10])

### Classify novel ORFs
orftable_annotated$floss_classification <- as.character(
    classifier_ref$classify(
        orftable_annotated$countRFP, orftable_annotated$floss
        )
        )

orftable_annotated <- orftable_annotated %>%
  mutate(floss_classification = case_when(
    floss_classification == "Good" ~ "Good",
     floss_classification == "Extreme" ~ "Extreme",
      is.na(floss_classification) ~ "ok"
  )
  )

orftable_annotated <- orftable_annotated %>%
  filter(floss_classification != "Extreme")

print(table(orftable_annotated$`Start codon`,orftable_annotated$`ORF type`))
## Outputs
### ORF identifications
if (!dir.exists("results/orf_identifications")) {
system("mkdir results/orf_identifications")
}

save(orftable_annotated, file="results/orf_identifications/ortable.rData")

### Fasta files
### Microproteins
microprotein_annotations <- orftable_annotated %>%
  filter(`ORF type` != "Annotated" ) %>%
  filter(`ORF type` == "Upstream" | `ORF type` == "Start overlap" |`ORF type` == "New") %>%
  filter(`Length (AAs)` < 100)

table(microprotein_annotations$`ORF type`,microprotein_annotations$`Start codon`)

cds_novel_orfs <- loadRegion(txdb_orf, part = "cds")

if (!dir.exists("results/orf_identifications")) {
system("mkdir results/fasta")
}

cgr_fasta <- FaFile("reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna")
orf_to_fasta(microprotein_annotations, cds_novel_orfs, "MicroProtein110823", cgr_fasta)

system("cat results/FASTA/UP000001075_10029.fasta results/FASTA/MicroProtein020823.fasta > results/FASTA/020823_UniProt_MicroProtein.fasta", intern = T)

## write a list of lncRNA encoded ORFs
## Only select the longest per transcript
selected_new_for_ts_de <- orftable_annotated %>%
filter(`ORF type` == "New") %>%
filter(str_detect(`ORF-RATER name`, "NR|XR")) %>%
group_by(`Transcript ID`) %>%
top_n(`Length (AAs)`, n=1)

if (!dir.exists("results/diff_trans_analysis")) {
system("mkdir results/diff_trans_analysis", intern = T,ignore.stderr = T, ignore.stdout = T)
}

write(selected_new_for_ts_de$`ORF-RATER name`, file="results/diff_trans_analysis/lncRNA_novel_orfs.txt")
