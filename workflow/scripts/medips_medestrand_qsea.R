 args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
sample_id  = args[1]
sample_bam = args[2]
bsgenome = args[3]
ispaired = as.logical(args[4])
medestrand_path = args[5]

# default setting
ws = 300

## loading libraries
library(MEDIPS)

## all chromosomes failed MEDIPS enrichment
## checking chr1-5 temporally

## loading corresponding genome and major chrs
## for testing dataset
if (bsgenome == "BSgenome.Scerevisiae.UCSC.sacCer3")
{
  library("BSgenome.Scerevisiae.UCSC.sacCer3")
  chr = "chrI"
  chr_enrich  = "chrI"
}

## hg19
if (bsgenome == "BSgenome.Hsapiens.UCSC.hg19")
{
  library("BSgenome.Hsapiens.UCSC.hg19")
  chr = c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  chr_enrich = "chr1"
}

## hg38
if (bsgenome == "BSgenome.Hsapiens.UCSC.hg38")
{
  library("BSgenome.Hsapiens.UCSC.hg38")

  chr = c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  chr_enrich = "chr1"

}


###################################
###         MEDIPS QC           ###
###################################

{
#############
## saturation
#############
saturation = MEDIPS.saturation(
  file = sample_bam,
  BSgenome = bsgenome,
  window_size = ws,
  paired = ispaired,
  chr.select = chr
)

## saturation plot
png(file = paste0("meth_qc_quant/", sample_id, "_saturation.png"), res = 300, width = 5, height = 5, units = "in")
MEDIPS.plotSaturation(saturationObj = saturation)
dev.off()


################
## CpG coverage
################
coverage = MEDIPS.seqCoverage(
  file = sample_bam,
  BSgenome = bsgenome,
  paired = ispaired,
  chr.select = chr
)

## coverage plot
png(file = paste0("meth_qc_quant/", sample_id, "_seqCoverage.png"), res = 300, width = 5, height = 5, units = "in")
MEDIPS.plotSeqCoverage(
  seqCoverageObj=coverage,
  type="pie",
  cov.level = c(0,1,2,3,4,5)
)
dev.off()

#################
## CpG enrichment
################
cpg_enrich = MEDIPS.CpGenrich(
  file = sample_bam,
  BSgenome = bsgenome,
  paired = ispaired,
  chr.select = chr_enrich
)

###################
## MEDIPS QC report
cov_sum = table(coverage$cov.res)
numberCpG = length(coverage$cov.res)

## fraction of CpG with more than 5 reads, not includes
idx_6 <- (names(cov_sum) != "0" & names(cov_sum) != "1" & names(cov_sum) != "2" &
          names(cov_sum) != "3" & names(cov_sum) != "4" & names(cov_sum) != "5")

## output qc matrix
medips_qc = data.frame(

  Sample = sample_id,

  ## saturation
  saturation.numberReads = saturation$numberReads,
  saturation.maxEstCor = saturation$maxEstCor[2],
  saturation.maxTruCor = saturation$maxTruCor[2],


  ## fraction of reads don't cover a CG
  coverage.fracReadsWoCpG = coverage$numberReadsWO / coverage$numberReads,

  ## fraction of CpG isn't covered by a read
  overage.fracCpGwoRead = cov_sum[names(cov_sum) == "0"] / numberCpG,

  ## fraction of CpG isn't covered by 1,2,3,4,5,>5 read(s)
  coverage.fracCpGw1Read = cov_sum[names(cov_sum) == "1"] / numberCpG,
  coverage.fracCpGw2Read = cov_sum[names(cov_sum) == "2"] / numberCpG,
  coverage.fracCpGw3Read = cov_sum[names(cov_sum) == "3"] / numberCpG,
  coverage.fracCpGw4Read = cov_sum[names(cov_sum) == "4"] / numberCpG,
  coverage.fracCpGw5Read = cov_sum[names(cov_sum) == "5"] / numberCpG,
  coverage.fracCpGgt5Reads = sum(cov_sum[idx_6]) / numberCpG,


  ## enrichment scores: closer to 1 indicating less enriched
  ## region
  enrichment.regions.CG = cpg_enrich$regions.CG,
  enrichment.regions.C = cpg_enrich$regions.C,
  enrichment.regions.G = cpg_enrich$regions.G,
  enrichment.regions.relH = cpg_enrich$regions.relH,
  enrichment.regions.GoGe = cpg_enrich$regions.GoGe,

  ## reference genome: will be consistants for same spciece
  enrichment.genome.CG = cpg_enrich$genome.CG,
  enrichment.genome.C = cpg_enrich$genome.C,
  enrichment.genome.G = cpg_enrich$genome.G,
  enrichment.genome.relH = cpg_enrich$genome.relH,
  enrichment.genome.GoGe = cpg_enrich$genome.GoGe,

  ## enrihment socres
  enrichment.relH = cpg_enrich$enrichment.score.relH,
  enrichment.GoGe = cpg_enrich$enrichment.score.GoGe

)

## consistent with multiQC sample X features output
write.table(medips_qc, file = paste0("meth_qc_quant/", sample_id, "_meth_qc.txt"),
            sep = "\t", row.names = F, quote = F)

}

###################################
###   MEDIPS quantification     ###
###################################
{
## MEDIPs set
mset = MEDIPS.createSet(
          file = sample_bam,
          BSgenome = bsgenome,
          paired = ispaired,
          window_size = ws,
          chr.select = chr
)

## coupling set
cset = MEDIPS.couplingVector(pattern = "CG", refObj = mset)

####################################################################
##  region CF(# of CpG), count, rpkm, rms:relative methylation level
meth =  MEDIPS.meth(
          MSet1 = mset,
          CSet = cset,
          CNV = F,
          MeDIP = T ,     ## T: rms
)

## rm redundant info for single sample
meth = meth[, 1:7]
colnames(meth)[c(4:7)] = c("CpG","count", "rpkm", "rms_medips")

########################################
## export coverage profile in wig format
# MEDIPS.exportWIG(Set = mset, CSet = cset, file = medips_rms_wig,
#                 format = "rms", descr = sample_id)

}



##################################################
## infer relative methylation level via MeDEStrand
#################################################
{
library(devtools)
devtools::load_all(medestrand_path)

## MeDEStrand set: reads mapped to the positive and negative DNA strands will be processed separately
medset <-  MeDEStrand.createSet(file = sample_bam,
            BSgenome = bsgenome,
            paired = ispaired,
            window_size = ws,
            chr.select = chr)

medcset <- MeDEStrand.countCG(pattern = 'CG', refObj = medset)

rms_medestrand <- MeDEStrand.binMethyl(
  MSetInput = medset,
  CSet = medcset,
  Granges = FALSE,
)

## adding to MEDIPS output
meth <- data.frame(meth, rms_medestrand)

}



##########################
### QSEA: under developing
##########################
{
library(qsea)

##################
## create qsea set
sample_info <- data.frame(sample_name = sample_id,
                          file_name = sample_bam,
                          group = sample_id)

qset = createQseaSet(sampleTable = sample_info,
                     BSgenome = bsgenome,
                     window_size = ws,
                     chr.select = chr)


#######################################
## count reads per genomic regions/bins
qset <- addCoverage(qset,
                    uniquePos = FALSE,     ## using dedup_bam
                    paired = ispaired)

## CNV: only fragments without CpG dinucleotides are considered for MEDIP .
qset = addCNV(qset,
              file_name="file_name",
              paired = ispaired,
              MeDIP = TRUE)

## add library factor and offset
qset <- addLibraryFactors(qset)

##################
## add CpG density
qset  <- addPatternDensity(qset, "CG", name="CpG")

## addOffset
qset <- addOffset(qset, "CpG")

###############################
## add enrichment parameters
## and plot enrichment profiles

## need to specify the windows for enrichment analysis
# Blind calibration: QSEA assumes that regions with low CpG density is 80% methylated
# on average, and regions within CpG islands are 25% methylated on average.
wd_idx <-  which(getRegions(qset)$CpG_density > 1 & getRegions(qset)$CpG_density < 10)
signal <- (15-getRegions(qset)$CpG_density[wd_idx])*.55/15+.25
signal <- matrix(signal, length(signal), 1)

qset <- addEnrichmentParameters(qset,
                                enrichmentPattern = "CpG",
                                windowIdx = wd_idx,
                                signal = signal)

################
## normal method
## nrpm: CNV normalized reads per million mappable reads
## beta: transformation to % methylation, posterior mean point estimator
## logitbeta: logit transformed beta values

nm = c("counts", "nrpm", "beta", "logitbeta")
qtb <- makeTable(qset,
                 samples = sample_id,
                 norm_methods = nm,
                 CNV = T)
colnames(qtb)[6:9] <- c("nrpm_qsea", "beta_qsea", "logitbeta_qsea", "CNV_qsea")


## QSEA will truncate the bin shorter than ws at the end chromation
## MEDIPS and MeDStrand will keep that bin by extending it to ws,
m_id <- paste(meth$chr, meth$start, meth$stop, sep = "_")
q_id <- paste(qtb$chr, qtb$window_start, qtb$window_end, sep = "_")
idx <- match(m_id, q_id)

meth_quant <- data.frame(meth[!is.na(idx), ], qtb[, 6:9])

## shared genomic regions and cpg count
#grange_cpg <- meth_quant[, 1:4]
#save(grange_cpg, file = paste0("meth_qc_quant/", sample_id, "_Granges_CpGs.Rdata"))

bin_id <- 1:nrow(meth_quant)
grange_cpg  <- data.frame(meth_quant[, 1:4], bin_id)
colnames(grange_cpg)[1] <- "#chr"
write.table(grange_cpg, file = paste0("meth_qc_quant/", sample_id, "_Granges_CpGs.bed"),
            row.names = F, col.names = T, sep = "\t", quote = F)


## export each quantification column as a separate file
L <- ncol(meth_quant)
for(i in 5:L) {
  out <- meth_quant[, i]
  feature <- colnames(meth_quant)[i]
  write.table(out, file = paste0("meth_qc_quant/", sample_id, "_", feature, ".txt"),
             row.names = F, col.names = sample_id, quote = F)
}

# print(object.size(meth_quant),units="auto")
save(meth_quant, file = paste0("meth_qc_quant/", sample_id, "_meth_quant.RData"))


}
