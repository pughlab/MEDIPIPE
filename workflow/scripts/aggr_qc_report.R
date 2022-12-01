args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
aggr_sample =     args[1]
ispaired = as.logical(args[2])
aggr_fastqc =     args[3]
aggr_sam_stats =  args[4]

if(ispaired){
aggr_frag_stats = args[5]
aggr_meth_qc =    args[6]
scr_dir =         args[7]
out_dir =         args[8]
} else {
aggr_meth_qc =    args[5]
scr_dir =         args[6]
out_dir =         args[7]
}

##
sample_aggr <- read.table(aggr_sample, header = T, sep = "\t")


###########################################################
## Reads QC (R1 only for both single- and paired-end reads)
###########################################################
{
fastqc <- read.table(aggr_fastqc, header = T, sep = "\t")

if(ispaired)
{
  raw_reads <- fastqc[seq(1, nrow(fastqc), 4), ]
  sample <- gsub("_R1", "", raw_reads$Sample)

  ## match with sample_aggr.tsv samples and group labels
  idx_g <- match(sample, sample_aggr$sample_id)
  sample <- sample[!is.na(idx_g)]
  group <- sample_aggr$group[idx_g[!is.na(idx_g)]]
  raw_reads <- raw_reads[!is.na(idx_g), ]

  prefilter_reads <- fastqc[seq(2, nrow(fastqc), 4), ]
  prefilter_reads$Sample <-  gsub("_R1_val_1", "", prefilter_reads$Sample)
  ## ensure the samples matched
  idx_m <- match(sample, prefilter_reads$Sample)
  prefilter_reads <- prefilter_reads[idx_m, ]

} else {
  #raw_reads <- fastqc[seq(1, nrow(fastqc), 2), ]
  ## only trimed fq reports aggreated for single-end
  raw_reads <- fastqc
  sample <- raw_reads$Sample

  ## match with sample_aggr.tsv samples and group labels
  idx_g <- match(sample, sample_aggr$sample_id)
  sample <- sample[!is.na(idx_g)]
  group <- sample_aggr$group[idx_g[!is.na(idx_g)]]
  raw_reads <- raw_reads[!is.na(idx_g), ]

  prefilter_reads <- raw_reads

  ##prefilter_reads$Sample <-  gsub("_val_1", "", prefilter_reads$Sample)
  ## ensure the samples matched
  ##idx_m <- match(sample, prefilter_reads$Sample)
  ##prefilter_reads <- fastqc[seq(2, nrow(fastqc), 2), ]
}

raw_reads_depth <- raw_reads$Total.Sequences
prefilter_reads_depth <- prefilter_reads$Total.Sequences
prefilter_reads_pct <- round(prefilter_reads_depth / raw_reads_depth, digits = 4)*100
prefilter_reads_GC <- prefilter_reads$X.GC
prefilter_reads_dedup_pct <- round(prefilter_reads$total_deduplicated_percentage, digits = 2)

qc_matrix <- data.frame(sample,  group, raw_reads_depth, prefilter_reads_depth,
                        prefilter_reads_pct, prefilter_reads_GC, prefilter_reads_dedup_pct)
}

###############
## alignment QC
###############
{
  sam_stats <- read.table(paste0(out_dir, "/", aggr_sam_stats), header = T, sep = "\t")

  raw_bam_stats <- sam_stats[seq(2, nrow(sam_stats), 2), ]
  raw_bam_stats$Sample <- gsub("_sorted", "", raw_bam_stats$Sample)

  ## match with qc_matrix order
  idx_m <- match(qc_matrix$sample, raw_bam_stats$Sample)
  raw_bam_stats <- raw_bam_stats[idx_m, ]

  ## after deduplication
  dedup_bam_stats <- sam_stats[seq(1, nrow(sam_stats), 2), ]
  dedup_bam_stats$Sample <- gsub("_dedup", "", dedup_bam_stats$Sample)

  # ensure sample matched
  idx_mm <- match(raw_bam_stats$Sample, dedup_bam_stats$Sample)
  dedup_bam_stats <-  dedup_bam_stats[idx_mm, ]

  ## against mapping reads
  mapped_reads_pct <- round(raw_bam_stats$reads_mapped_percent, digits = 2)

  ## against all raw read
  if(ispaired)
  {
    ## counting read pairs
    usable_reads_depth <- dedup_bam_stats$sequences/2
    usable_reads_depth_pct <- round(usable_reads_depth  / qc_matrix$raw_reads_depth, digits = 4)*100
  } else {
    ## counting single-end reads
    usable_reads_depth <- dedup_bam_stats$sequences
    usable_reads_depth_pct <- round(usable_reads_depth  / qc_matrix$raw_reads_depth, digits = 4)*100
  }

  qc_matrix <- data.frame(qc_matrix, mapped_reads_pct, usable_reads_depth, usable_reads_depth_pct)



}

###############
## Fragment QC
###############
if(ispaired){
  frag_stats <- read.table(paste0(out_dir, "/", aggr_frag_stats), header = T, sep = "\t")
  frag_stats$Sample <- gsub("_dedup_FR", "", frag_stats$Sample)

  ## match with qc_matrix order
  idx_m <- match(qc_matrix$sample, frag_stats$Sample)
  frag_stats <- frag_stats[idx_m, ]

  fragment_size_mode <- frag_stats$MODE_INSERT_SIZE
  fragment_size_mean <- round(frag_stats$MEAN_INSERT_SIZE, digits = 0)
  fragment_size_median <- frag_stats$MEDIAN_INSERT_SIZE

  ## WIDTH_OF_80_PERCENT
  ## The "width" of the bins, centered around the median, that encompass 80% of all read pairs.
  fragment_size_for_80pct_reads <- frag_stats$MEDIAN_INSERT_SIZE + frag_stats$WIDTH_OF_80_PERCENT/2


  qc_matrix <- data.frame(qc_matrix, fragment_size_mode,
                          fragment_size_mean, fragment_size_median,
                          fragment_size_for_80pct_reads)

  }

###########
## MeDIP QC
###########
{
  meth_qc <- read.table(paste0(out_dir, "/", aggr_meth_qc), header = T, sep = "\t")
  idx_m <- match(qc_matrix$sample, meth_qc$Sample)
  meth_qc <- meth_qc[idx_m, ]

  saturation_maxEstCor <- round(meth_qc$saturation.maxEstCor, digits = 2)
  coverage_pctReadsWoCpG <- round(meth_qc$coverage.fracReadsWoCpG, digits = 4) * 100
  coverage_pctCpGwoRead <- round(meth_qc$coverage.fracCpGwoRead, digits = 4) * 100
  coverage_pctCpGwo1Read <- round(meth_qc$coverage.fracCpGw1Read, digits = 4) * 100
  coverage_pctCpGgt5Read <- round(meth_qc$coverage.fracCpGgt5Reads, digits = 4) * 100
  enrichment_relH <- round(meth_qc$enrichment.relH, digits = 2)
  enrichment_GoGe <- round(meth_qc$enrichment.GoGe, digits = 2)

  qc_matrix <- data.frame(qc_matrix, saturation_maxEstCor, coverage_pctReadsWoCpG,
                          coverage_pctCpGwoRead, coverage_pctCpGwo1Read, coverage_pctCpGgt5Read,
                          enrichment_relH, enrichment_GoGe)
}

##########################
## Generate HTML QC report
##########################
{

write.csv(qc_matrix, "aggregated/aggr_qc_report.csv", row.names = F)

## metrics correlatioin
library(corrplot)
mx <- as.matrix(qc_matrix[, 3:ncol(qc_matrix)])
M <- cor(mx)
pdf(file = "aggregated/aggr_qc_corrplot.pdf", width = 8, height = 8)
corrplot(M)
dev.off()

library(rmarkdown)
render(paste0(scr_dir, "/workflow/scripts/aggr_qc_report.Rmd"), output_dir = "aggregated",
       params = list(readin = paste0(out_dir, "/aggregated/aggr_qc_report.csv"), ispaired = ispaired))
}
