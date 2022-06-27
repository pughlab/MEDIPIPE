## This script was developed based on : 
## https://github.com/cancer-genomics/delfi_scripts
args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
sample_id  = args[1]
sample_bam = args[2]      ## should be indexed paired-end bam file
bsgenome = args[3]

## for testing 
if(FALSE){
rm(list = ls())
setwd("/Users/yong/OneDrive - UHN/Projects/TCGE/Hansen_group/Fragment_ratio")
sample_id = "test"
sample_bam = "test.bam"          
bsgenome = "BSgenome.Hsapiens.UCSC.hg19"
}

if (bsgenome == "BSgenome.Hsapiens.UCSC.hg19")
{
  library("BSgenome.Hsapiens.UCSC.hg19")
  Hsapiens <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  ## load gaps and filters
  load("gaps_filters_hg19.rdata")
  ## load AB: hi_AB_compartments
  load("AB_hg19.rdata")
  
} else if (bsgenome == "BSgenome.Hsapiens.UCSC.hg38") {
  library("BSgenome.Hsapiens.UCSC.hg38")
  Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  ## load gaps and filters
  load("gaps_filters_hg38.rdata")
  ## load AB: hi_AB_compartments
  load("AB_hg38.rdata")
}

## required packages 
library(GenomicAlignments)
library(biovizBase)            ##GCcontent function 
library(dplyr)

#########################
## GC correction function
gc.correct <- function(coverage, bias) {
    i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
    coverage.trend <- loess(coverage ~ bias)
    coverage.model <- loess(predict(coverage.trend, i) ~ i)
    coverage.pred <- predict(coverage.model, bias)
    coverage.corrected <- coverage - coverage.pred + median(coverage)
}

####################
## Fragment from bam
# parameters  using for function of readGAlignmentPairs
param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE,
                                         isUnmappedQuery = FALSE), mapqFilter = 30)
galp <- readGAlignmentPairs(sample_bam, param = param)

## autosome only
fragments <- granges(keepSeqlevels(galp, paste0("chr", 1:22), pruning.mode="coarse"),
                     on.discordant.seqnames="drop")

## filter out extreme values 
w.all <- width(fragments)
q.all <- quantile(w.all, c(0.001, 0.999))
fragments <- fragments[which(w.all > q.all[1] & w.all < q.all[2])]

## add gc content 
fragments$gc <- GCcontent(Hsapiens, unstrand(fragments))

######################################
## remove fragments in filters regions
fragments <- fragments[-queryHits(findOverlaps(fragments, filters))];

## focusing on fragment length in [100, 220]
w.all <- width(fragments)
fragments <- fragments[which(w.all >= 100 & w.all <= 220)]
w <- width(fragments)

## group fragments by their length
frag.list <- split(fragments, w)


################################################################################
## HiC_AB_Compartments
## data source (hg19): https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt
#library(RCurl)
#ABurl <- getURL('https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
# liftover to hg38
#AB <- read.table("hic_compartments_100kb_ebv_2014_liftover_hg38.txt", header=TRUE)
#AB <- makeGRangesFromDataFrame(AB, keep.extra.columns=TRUE)
#save(AB, file = "AB_hg38.rdata")

## filter out gaps 
chromosomes <- GRanges(paste0("chr", 1:22), IRanges(0, seqlengths(Hsapiens)[1:22]))
tcmeres <- gaps[grepl("centromere|telomere", gaps$type)]

## remove centromere|telomere
arms <- GenomicRanges::setdiff(chromosomes, tcmeres)

arms <- arms[-c(25,27,29,41,43)]   ## ??
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
arms$arm <- armlevels

## remove AB bins [100k] overlapped with Gaps
AB <- AB[-queryHits(findOverlaps(AB, gaps))]
AB <- AB[queryHits(findOverlaps(AB, arms))]       ## reomved 5 arms
AB$arm <- armlevels[subjectHits(findOverlaps(AB, arms))]

seqinfo(AB) <- seqinfo(Hsapiens)[seqlevels(seqinfo(AB))]
AB <- trim(AB)
AB$gc <- GCcontent(Hsapiens, AB)


######### 
olaps <- findOverlaps(fragments, AB)
bin.list <- split(fragments[queryHits(olaps)], subjectHits(olaps))

## mean GC of the fragments overlapped with specific bin
#bingc <- rep(NA, length(bin.list))
bingc <- rep(NA, length(AB))   ## preserve AB bins without coverage, ensure same AB dim
bingc[unique(subjectHits(olaps))] <- sapply(bin.list, function(x) mean(x$gc))


## count the number of fragments with certain length in each AB bins
counts <- sapply(frag.list, function(x) countOverlaps(AB, x))

## make sure the counts colnames range [100, 220]
if(min(w) > 100) {
  m0 <- matrix(0, ncol=min(w) - 100, nrow=nrow(counts),
               dimnames=list(rownames(counts), 100:(min(w)-1)))
  counts <- cbind(m0, counts)
}

#########
## output
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## stats of fragments after filtering
modes <- Mode(w)
medians <- median(w)
q25 <- quantile(w, 0.25)
q75 <- quantile(w, 0.75)

## cnt and ratio of short [100-150] and long [151: 220] fragments 
short <- rowSums(counts[,1:51])
long <- rowSums(counts[,52:121])
ratio <- short/long
short.corrected=gc.correct(short, bingc);
long.corrected=gc.correct(long, bingc)
nfrags.corrected=gc.correct(short+long, bingc)

ratio.corrected <- rep(NA, length(ratio));
ratio.corrected[is.finite(ratio)]=gc.correct(ratio[is.finite(ratio)], bingc[is.finite(ratio)]);
AB$short <- short
AB$long <- long
AB$ratio <- short/long
AB$nfrags <- short+long
AB$short.corrected <- short.corrected
AB$long.corrected <- long.corrected
AB$nfrags.corrected <- nfrags.corrected
AB$ratio.corrected <- ratio.corrected

AB$mode <- modes
AB$mean <- round(mean(w), 2)
AB$median <- medians
AB$quantile.25 <- q25
AB$quantile.75 <- q75
AB$frag.gc <- bingc

## combine AB and counts matrix
for(i in 1:ncol(counts)) elementMetadata(AB)[,colnames(counts)[i]] <- counts[,i]
saveRDS(AB, paste0(sample_id, "_100kb_fragment_profile.RDS"))

## transfer to data.frame
dat <- as_tibble(AB)
dat$arm <- factor(dat$arm, levels=armlevels)

for (bin_size in c(10, 50))   ## bin size to be 1mb (10*100kb) and 5mb 
{
  ## combine adjacent 100kb bins to form bin_size bins. We count starting from
  ## the telomeric end and remove the bin closest to the centromere if it is
  ## smaller than 5mb.
  dat_s <-  dat %>% 
            group_by(arm) %>%
            mutate(combine = ifelse(grepl("p", arm), ceiling((1:length(arm))/bin_size),
                             ceiling(rev((1:length(arm))/bin_size))))
  
  ## summarize by new bin_size 
  dat_ss <- dat_s %>% 
            group_by(seqnames, arm, combine) %>%
            summarize(short2=sum(short),
              long2=sum(long),
              short.corrected2=sum(short.corrected),
              long.corrected2=sum(long.corrected),
              hic.eigen=mean(eigen),
              gc=mean(C.G),
              ratio2=mean(ratio),
              ratio.corrected2=mean(ratio.corrected),
              nfrags2=sum(nfrags),
              nfrags.corrected2=sum(nfrags.corrected),
              domain = median(as.integer(domain)),
              short.var=var(short.corrected),
              long.var=var(long.corrected),
              nfrags.var=var(nfrags.corrected),
              mode_size=unique(mode),
              mean_size=unique(mean),
              median_size=unique(median),
              q25_size=unique(quantile.25),
              q75_size=unique(quantile.75),
              start=start[1],
              end=rev(end)[1],
              binsize = n())
  ## filter out bins less than selected bin_size
  dat_ss <- dat_ss %>% filter(binsize==bin_size);
  saveRDS( dat_ss, paste0(sample_id, "_", bin_size, "_100kb_fragment_profile.RDS"))
  
  ## print out corrected frament ratio 
  pr_out <- data.frame(dat_ss$seqnames, dat_ss$start, dat_ss$end, dat_ss$ratio.corrected2) 
  write.table(pr_out, file = paste0(sample_id, "_", bin_size, "_100kb_fragment_profile_GC_corrected_Ratio.txt"),
              quote = F, row.names = F, col.names = F)
  
}


