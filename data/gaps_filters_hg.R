## get human genome Gaps and ENOCDE black list regions
## the code modified by Yong basd on:
## https://github.com/cancer-genomics/delfi_scripts/blob/master/00-filtered_regions.r

library(rtracklayer)
library(tidyverse)

## For hg19
if(FALSE){
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
## Blacklist source: https://storage.googleapis.com/encode-pipeline-genome-data/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz
blacklist <- "wgEncodeDacMapabilityConsensusExcludable.bed"
genome <- "hg19"
}

## for hg38
if(TRUE){
library(BSgenome.Hsapiens.UCSC.hg38)
Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
genome <- "hg38"
## Blacklist source: https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
blacklist <- "ENCFF356LFX.bed"
}

#######
## gaps
mySession <- browserSession()
genome(mySession) <- genome
gaps <- getTable(ucscTableQuery(mySession, track="gap"))
gaps <- GRanges(gaps$chrom, IRanges(gaps$chromStart,
                     gaps$chromEnd),type=gaps$type)

gaps <- keepSeqlevels(gaps, paste0("chr", c(1:22, "X", "Y")),
                           pruning.mode="coarse")
seqinfo(gaps) <- seqinfo(Hsapiens)[seqlevels(gaps), ]

#### encdoe black list 
blacklisted.tib <- read_tsv(blacklist, col_names=c("seqnames", "start","end"))
blacklisted.tib <- blacklisted.tib %>% mutate(start=start+1)
filters <- makeGRangesFromDataFrame(blacklisted.tib, keep.extra.columns=TRUE)
filters <- keepSeqlevels(filters, paste0("chr", c(1:22, "X", "Y")), pruning.mode="coarse")
seqinfo(filters) <- seqinfo(Hsapiens)[seqlevels(filters), ]

save(gaps, filters, file = paste0("gaps_filters_", genome, ".rdata"))
     
     