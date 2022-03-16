args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
sampleid = args[1]
bam_file = args[2]
bsgenome = args[3]
ispaired = as.logical(args[4])
satu_png = args[5]
covr_png = args[6]
qc_repor = args[7]
meth_qua = args[8]
rpkm_wig = args[9]
medestrand = args[11]

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
  chr = c(paste0("chr", 1:22), "chrX", "chrY")
  chr_enrich = "chr1"
}

## hg38
if (bsgenome == "BSgenome.Hsapiens.UCSC.hg38")
{
  library("BSgenome.Hsapiens.UCSC.hg38")

  chr = c(paste0("chr", 1:22), "chrX", "chrY")
  chr_enrich = "chr1"

}


###################################
###         MEDIPS QC           ###
###################################

{
#############
## saturation
saturation = MEDIPS.saturation(
  file = bam_file,
  BSgenome = bsgenome,
  window_size = ws,
  paired = ispaired,
  chr.select = chr
)

## saturation plot
png(satu_png, res = 300, width = 5, height = 5, units = "in")
MEDIPS.plotSaturation(saturationObj = saturation)
dev.off()



################
## CpG coverage
################

coverage = MEDIPS.seqCoverage(
  file = bam_file,
  BSgenome = bsgenome,
  paired = ispaired,
  chr.select = chr
)

## coverage plot
png(file = covr_png, res = 300, width = 5, height = 5, units = "in")
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
  file = bam_file,
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

  Sample = sampleid,

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

write.table(medips_qc, file = qc_repor,
            sep = "\t", row.names = F, quote = F)

}

###################################
###   MEDIPS quantification     ###
###################################
{
## MEDIPs set
mset = MEDIPS.createSet(
          file = bam_file,
          BSgenome = bsgenome,
          paired = ispaired,
          window_size = ws,
          chr.select = chr
)

#############################################################
## methylation profile: region CF(# of CpG), count, rpkm, rms
if (FALSE) {
## coupling set
cset = MEDIPS.couplingVector(pattern = "CG", refObj = mset)

meth =  MEDIPS.meth(
          MSet1 = mset,
          CSet = cset,
          CNV = F,
          MeDIP = T ,           ## calculate rms: relative methylation levels
)

## rm redundant info for single sample
meth = meth[, 1:7]
colnames(meth)[5:7] = c("count", "rpkm", "rms")
write.table(meth, file = meth_qua,
            row.names = F, quote = F, sep ="\t")
}

##########################
## without calculating rms
meth =  MEDIPS.meth(
          MSet1 = mset,
          CNV = F,
          MeDIP = F ,
)

## rm redundant info for single sample
meth = meth[, 1:5]
colnames(meth)[4:5] = c("count", "rpkm")
write.table(meth, file = meth_qua,
            row.names = F, quote = F, sep ="\t")

########################################
## export coverage profile in wig format
MEDIPS.exportWIG(Set = mset, file = rpkm_wig,
                 format = "rpkm", descr = sampleid)

}



############################################
## infer abosulate m6A leves via MeDEStrand
###########################################
{
library(devtools)
devtools::load_all(medestrand)


## MeDEStrand set: reads mapped to the positive and negative DNA strand are processed separately
medset <-  MeDEStrand.createSet(file = bam_file,
            BSgenome = bsgenome,
            paired = ispaired,
            window_size = ws,
            chr.select = chr)

medcset <- MeDEStrand.countCG(pattern = 'CG', refObj = medset)

meth_abs <- MeDEStrand.binMethyl(
  MSetInput = medset,
  CSet = medcset,
  Granges = FALSE,
)

saveRDS(meth_abs, file = paste0("meth_quant/", sampleid, "_meth_abs.RDS"))
}

##########################
### QSEA: under developing
##########################
if(FALSE)
{
library(qsea)

##################
## create qsea set
sample_info <- data.frame(sample_name = sample_id,
                          file_name = bam,
                          group = sample_id)

qs = createQseaSet(sampleTable = sample_info,
                   BSgenome = bsgenome,
                   window_size = ws,
                   chr.select = chr)


#######################################
## count reads per genomic regions/bins
qs <- addCoverage(qs,
            uniquePos = FALSE,   ## removed PCR duplicates in previous step
            paired = ispaired)


##################
## add CpG density
qs  <- addPatternDensity(qs, "CG", name="CpG")

#################################
## add library factor and offset
qs <- addLibraryFactors(qs)
qs <- addOffset(qs, "CpG", maxPatternDensity=0.7)


###############################
## add enrichment parameters
## and plot enrichment profiles

## need to specify the windows for enrichment analysis
# QSEA assumes that regions with low CpG density is 80% methylated
# on average, and regions within CpG islands are 25% methylated on average.

wd_idx <-  which(getRegions(qs)$CpG_density > 1 & getRegions(qs)$CpG_density < 10)
signal <- (15-getRegions(qs)$CpG_density[wd_idx])*.55/15+.25
signal <- matrix(signal, length(signal), 1)

qs <- addEnrichmentParameters(qs,
                              enrichmentPattern = "CpG",
                              windowIdx = wd_idx,
                              signal = signal)

}
