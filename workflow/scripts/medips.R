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

# default setting
ws = 300

## loading libraries
library(MEDIPS)

## loading corresponding genome and major chrs
## for testing dataset
if (bsgenome == "BSgenome.Scerevisiae.UCSC.sacCer3")
{
  library("BSgenome.Scerevisiae.UCSC.sacCer3")
  chr = "chrI"
}

## hg19
if (bsgenome == "BSgenome.Hsapiens.UCSC.hg19")
{
  library("BSgenome.Hsapiens.UCSC.hg19")
  chr = c(pate0("chr", 1:22), "chrX", "chrY")
}

## hg38
if (bsgenome == "BSgenome.Hsapiens.UCSC.hg38")
{
  library("BSgenome.Hsapiens.UCSC.hg38")
  chr = c(pate0("chr", 1:22), "chrX", "chrY")
}


###################################
###         MEDIPS QC           ###
###################################

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
cpg_enrich = MEDIPS.CpGenrich(
  file = bam_file,
  BSgenome = bsgenome,
  paired = ispaired,
  chr.select = chr
)

###################
## MEDIPS QC report
cov_sum = table(coverage$cov.res)
numberCpG = length(coverage$cov.res)

medips_qc = data.frame(

  sample = sampleid,

  saturation.numberReads = saturation$numberReads,
  saturation.maxEstCor = saturation$maxEstCor[2],
  saturation.maxTruCor = saturation$maxTruCor[2],

  ## fraction of reads don't cover a CG
  coverage.fracReadsWoCpG = coverage$numberReadsWO / coverage$numberReads,
  ## fraction of CpG isn't covered by a read
  coverage.fracCpGwoRead = cov_sum[1] / numberCpG,
  ## fraction of CpG with more than 5 reads, not includes 5
  coverage.fracCpGgt5Reads = sum(cov_sum[7:length(cov_sum)]) / numberCpG,

  ## enrichment scores: closer to 1 indicating less enriched
  enrichment.relH = cpg_enrich$enrichment.score.relH,
  enrichment.GoGe = cpg_enrich$enrichment.score.GoGe
)

write.table(medips_qc, file = qc_repor,
            sep = "\t", row.names = F, quote = F)

###################################
###   MEDIPS quantification     ###
###################################

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
