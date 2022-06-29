## get human genome Gaps and ENOCDE black list regions
## the code modified by Yong basd on:
## https://github.com/cancer-genomics/delfi_scripts/blob/master/00-filtered_regions.r

library(rtracklayer)
library(tidyverse)

## For hg19
{
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
## Blacklist source: https://storage.googleapis.com/encode-pipeline-genome-data/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz
blacklist <- "wgEncodeDacMapabilityConsensusExcludable.bed"
genome <- "hg19"

#######
## gaps
mySession <- browserSession()
genome(mySession) <- genome

## there are no centromeres for hg38 gaps!!!
gaps <- getTable(ucscTableQuery(mySession, table="gap"))
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
}

###########
## for hg38
{
library(BSgenome.Hsapiens.UCSC.hg38)
Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
genome <- "hg38"
## Blacklist source: https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
blacklist <- "ENCFF356LFX.bed"

#######
## gaps
mySession <- browserSession()
genome(mySession) <- genome

gaps <- getTable(ucscTableQuery(mySession, table="gap"))
gaps <- GRanges(gaps$chrom, IRanges(gaps$chromStart,
                                    gaps$chromEnd),type=gaps$type)

## there are no centromeres for hg38 gaps!!!
## download centromeres from:
## http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1390519713_pXF2zmUJbgVylZHcrke8Mu3VHvax&clade=mammal&org=Human&db=hg38&hgta_group=map&hgta_track=centromeres&hgta_table=0&hgta_regionType=genome&position=chrX%3A15%2C560%2C138-15%2C602%2C945&hgta_outputType=primaryTable&hgta_outFileName=
## merge to single centromere per chromosome by taking the [min, max]
cent_all <- read.table("hg38_centromeres.bed", header = F)
chrs <- unique(cent_all$V1)
cent_m <- data.frame()
for(i in 1:length(chrs))
{
  cent_m[i, 1] <- chrs[i]
  idx <- cent_all$V1 == chrs[i]
  cent_m[i, 2] <- min(cent_all$V2[idx])
  cent_m[i, 3] <- max(cent_all$V3[idx])
}
cent_mg <- GRanges(cent_m$V1, IRanges(cent_m$V2, cent_m$V3), type="centromere")

gaps <- c(gaps, cent_mg)
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
}
