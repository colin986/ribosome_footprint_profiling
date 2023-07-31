# functions for ORF filtering

count_transcript_alignments <- function(bam_data) {
  transcript_granges <- GRanges(
    seqnames = names(seqlengths(bam_data)),
    ranges = IRanges(1, end = seqlengths(bam_data)),
    strand = "*")
  
  transcript_counts <- countOverlaps(transcript_granges, bam_data)
  
  counts <- data.frame(transcript = as.character(transcript_granges@seqnames),
                       counts = transcript_counts)
  
  return(counts)
}

# this function was taken from ingolia et al. Floss paper
library(zoo)
flossCutoffs <- function(nreads, score) {
  rollwidth <- 200
  
  unsorted <- data.frame(nreads=nreads, score=score)  
  fr <- unsorted[order(unsorted$nreads),]
  nr <- rollapply(fr$nreads, width=rollwidth, function(xs) { mean(xs, na.rm=T) })
  medsc <- rollapply(fr$score, width=rollwidth, function(xs) { median(xs, na.rm=T) })
  q1 <- rollapply(fr$score, width=rollwidth, function(xs) { quantile(xs, probs=0.25, na.rm=T) })
  q3 <- rollapply(fr$score, width=rollwidth, function(xs) { quantile(xs, probs=0.75, na.rm=T) })
  rawExtreme <- q3 + 3.0*(q3-q1)
  extLoess <- loess(extreme ~ nread, data.frame(nread=log10(nr), extreme=rawExtreme))
  
  classify <- function(nr, sc) {
    extCutoff <- predict(extLoess, log10(nr))
    factor(
      if (is.na(extCutoff)) {
        NA
      } else if (sc > extCutoff) {
        "Extreme"      
      } else {
        "Good"
      },
      levels=c("Good", "Extreme")
    )
  }
  
  list( l10nreads = log10(nr),
        scoreMedian = medsc, scoreQ3 = q3, scoreExtreme = rawExtreme,
        extremeLoess = extLoess,
        isExtreme = function(nr, sc) { sc > predict(extLoess, log10(nr)) },
        classify = Vectorize(classify)
  )
}

orf_to_fasta <- function(orftable, cds_granges, name, fasta) {
  
  standard_code_with_alt_init_codon <- GENETIC_CODE
  attr(standard_code_with_alt_init_codon, "alt_init_codons") <- c("CTG", "GTG","TTG")
  
  orf_granges <- cds_granges[names(cds_granges) %in% orftable$`ORF-RATER name`]
  
  orf_transcript_seqs <- extractTranscriptSeqs(fasta, 
                                               orf_granges)
  
  orf_protein_seqs <- Biostrings::translate(orf_transcript_seqs,genetic.code = standard_code_with_alt_init_codon)
  
  
  output_path = paste0("../../results/FASTA/",name,".fasta")
  
  writeXStringSet(orf_protein_seqs, filepath = output_path)
}


overlapping <- function(x){
  if (length(unique(x$`Transcript stop position`)) < length(x$`Transcript stop position`)) {
    n_occur <- data.frame(table(x$`Transcript stop position`))
    x_d <- x[duplicated(x$`Transcript stop position`),] 
    x_d <- x_d[which.max(x_d$`Transcript stop position`-x_d$`Transcript start position`),] 
    if (length(n_occur$Freq == 1) > 0) { 
      nd_tstop <- n_occur[n_occur$Freq <2,1] 
      x_nd <- x[x$`Transcript stop position` %in% nd_tstop,] 
      x_final <- rbind(x_nd, x_d) } 
    else { 
      x_final <- x_d
    }
  } else  {
    x_final <- x
  }
  x_final
}
