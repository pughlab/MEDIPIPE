args = commandArgs(trailingOnly=TRUE)

## assigning commandArgs
sample_id  = args[1]
sample_bam = args[2]
ispaired = as.logical(args[3])
bsgenome = args[4]
bsgenome_pkg = args[5]

# default setting
ws = 300

## loading libraries
library(MEDIPS)

## all chromosomes failed MEDIPS enrichment
## checking chr1 only temporally, which will lead to slightly higher enrichment scores
## Problem has beed fixed by the function MEDIPS.CpGenrichNew below

## loading corresponding genome and major chrs
## for spike-in

if (is.element(bsgenome, installed.packages()[,1]))
{
  library("BSgenome.Athaliana.BAC.F19K16.F24B22")
  chr = c("AC011717.6", "AL132957.1")
} else {
  ## install custom spike-in BSgenome pkg
  install.packages(bsgenome_pkg, repos = NULL, type="source")

  library("BSgenome.Athaliana.BAC.F19K16.F24B22")
  chr = c("AC011717.6", "AL132957.1")
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
  uniq = 0,
  chr.select = chr
)

## saturation plot
png(file = paste0("meth_qc_quant_spikein/", sample_id, "_saturation.png"), res = 300, width = 5, height = 5, units = "in")
MEDIPS.plotSaturation(saturationObj = saturation)
dev.off()


################
## CpG coverage
################
coverage = MEDIPS.seqCoverage(
  file = sample_bam,
  BSgenome = bsgenome,
  paired = ispaired,
  uniq = 0,
  chr.select = chr
)

## coverage plot
png(file = paste0("meth_qc_quant_spikein/", sample_id, "_seqCoverage.png"), res = 300, width = 5, height = 5, units = "in")
MEDIPS.plotSeqCoverage(
  seqCoverageObj=coverage,
  type="pie",
  cov.level = c(0,1,2,3,4,5)
)
dev.off()

#################
## CpG enrichment
#################

#########################################################
## rewritten enrichment function based on SPPearce's code
## https://github.com/chavez-lab/MEDIPS/issues/2b
MEDIPS.CpGenrichNew <-function(file=NULL, BSgenome=NULL, extend=0, shift=0, uniq=1e-3, chr.select=NULL, paired=F){

	## Proof correctness....
  if(is.null(BSgenome)){stop("Must specify a BSgenome library.")}

  ## Read region file
  fileName=unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
  path=paste(unlist(strsplit(file, "/"))[1:(length(unlist(strsplit(file, "/"))))-1], collapse="/")
  if(path==""){path=getwd()}
  if(!fileName%in%dir(path)){stop(paste("File", fileName, " not found in", path, sep =" "))}

  #dataset = get(ls(paste("package:", BSgenome, sep = "")))
  dataset = get(ls(paste("package:", BSgenome, sep = ""))[1])   ## first element

  if(!paired){GRange.Reads = getGRange(fileName, path, extend, shift, chr.select, dataset, uniq)}	else{GRange.Reads = getPairedGRange(fileName, path, extend, shift, chr.select, dataset, uniq)}

  ## Sort chromosomes
  library(gtools)             ## need to load gtools
  if(length(unique(seqlevels(GRange.Reads)))>1){chromosomes=mixedsort(unique(seqlevels(GRange.Reads)))}
  if(length(unique(seqlevels(GRange.Reads)))==1){chromosomes=unique(seqlevels(GRange.Reads))}

  ## Get chromosome lengths for all chromosomes within data set.
  cat(paste("Loading chromosome lengths for ",BSgenome, "...\n", sep=""))

  chr_lengths=as.numeric(seqlengths(dataset)[chromosomes])

  ranges(GRange.Reads) <- restrict(ranges(GRange.Reads),+1)

  ##Calculate CpG density for regions
  total=length(chromosomes)
  cat("Calculating CpG density for given regions...\n")

  readsChars <- unlist(getSeq(dataset, GRange.Reads, as.character=TRUE))

  regions.CG = sum(vcountPattern("CG",readsChars))
  regions.C  = sum(vcountPattern("C",readsChars))
  regions.G  = sum(vcountPattern("G",readsChars))
  all.genomic= sum(width(readsChars))

  nReads <- length(readsChars)

  regions.relH=as.numeric(regions.CG)/as.numeric(all.genomic)*100
  regions.GoGe=(as.numeric(regions.CG)*as.numeric(all.genomic))/(as.numeric(regions.C)*as.numeric(regions.G))

 CG <- DNAStringSet("CG")
  pdict0 <- PDict(CG)
  params <- new("BSParams", X = dataset, FUN = countPDict, simplify = TRUE, exclude = c("rand", "chrUn"))
  genome.CG=sum(bsapply(params, pdict = pdict0))
  params <- new("BSParams", X = dataset, FUN = alphabetFrequency, exclude = c("rand", "chrUn"), simplify=TRUE)
  alphabet=bsapply(params)
  genome.l=sum(as.numeric(alphabet))
  genome.C=as.numeric(sum(alphabet[2,]))
  genome.G=as.numeric(sum(alphabet[3,]))
  genome.relH=genome.CG/genome.l*100
  genome.GoGe=(genome.CG*genome.l)/(genome.C*genome.G);

  ##Calculate CpG density for reference genome

  enrichment.score.relH=regions.relH/genome.relH
  enrichment.score.GoGe=regions.GoGe/genome.GoGe

  gc()
  return(list(genome=BSgenome, regions.CG=regions.CG, regions.C=regions.C, regions.G=regions.G, regions.relH=regions.relH, regions.GoGe=regions.GoGe, genome.C=genome.C, genome.G=genome.G, genome.CG=genome.CG, genome.relH=genome.relH, genome.GoGe=genome.GoGe, enrichment.score.relH=enrichment.score.relH, enrichment.score.GoGe=enrichment.score.GoGe))
}

## apply new enrichment function
cpg_enrich = MEDIPS.CpGenrichNew(
  file = sample_bam,
  BSgenome = bsgenome,
  paired = ispaired,
  uniq = 0,
  chr.select = chr
)


###################
## MEDIPS QC report
cov_sum = table(coverage$cov.res)
numberCpG = length(coverage$cov.res)

## fraction of CpG with more than 5 reads, not includes
idx_6 <- (names(cov_sum) != "0" & names(cov_sum) != "1" & names(cov_sum) != "2" &
          names(cov_sum) != "3" & names(cov_sum) != "4" & names(cov_sum) != "5")

## CpGs could be not covered by certain number of reads for spikeins
## fraction of CpG isn't covered by a read
if(length(cov_sum[names(cov_sum) == "0"]) == 0){
coverage.fracCpGwoRead = 0} else {
coverage.fracCpGwoRead = cov_sum[names(cov_sum) == "0"] / numberCpG}

## fraction of CpG isn't covered by 1,2,3,4,5,>5 read(s)
if(length(cov_sum[names(cov_sum) == "1"]) == 0){
coverage.fracCpGw1Read = 0} else {
coverage.fracCpGw1Read = cov_sum[names(cov_sum) == "1"] / numberCpG}

if(length(cov_sum[names(cov_sum) == "2"]) == 0){
coverage.fracCpGw2Read = 0} else {
coverage.fracCpGw2Read = cov_sum[names(cov_sum) == "2"] / numberCpG}

if(length(cov_sum[names(cov_sum) == "3"]) == 0){
coverage.fracCpGw3Read = 0} else {
coverage.fracCpGw3Read = cov_sum[names(cov_sum) == "3"] / numberCpG}

if(length(cov_sum[names(cov_sum) == "4"]) == 0){
coverage.fracCpGw4Read = 0} else {
coverage.fracCpGw4Read = cov_sum[names(cov_sum) == "4"] / numberCpG}

if(length(cov_sum[names(cov_sum) == "5"]) == 0){
coverage.fracCpGw5Read = 0} else {
coverage.fracCpGw5Read = cov_sum[names(cov_sum) == "5"] / numberCpG}


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
  coverage.fracCpGwoRead,

  ## fraction of CpG isn't covered by 1,2,3,4,5,>5 read(s)
  coverage.fracCpGw1Read,
  coverage.fracCpGw2Read,
  coverage.fracCpGw3Read,
  coverage.fracCpGw4Read,
  coverage.fracCpGw5Read,
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
write.table(medips_qc, file = paste0("meth_qc_quant_spikein/", sample_id, "_meth_qc.txt"),
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
	  uniq = 0,
          chr.select = chr
)

## coupling set
cset = MEDIPS.couplingVector(pattern = "CG", refObj = mset)
CpG = cset@genome_CF
####################################################################
##  region CF(# of CpG), count, rpkm, rms:relative methylation level
meth =  MEDIPS.meth(
          MSet1 = mset,
          CSet = cset,
          CNV = F,
          MeDIP = F ,     ## Failed rms estimation for spikeins
)

## rm redundant info for single sample
meth = data.frame(meth[, 1:5], CpG)
colnames(meth)[c(4:5)] = c("count", "rpkm")

## Granges and CpG counts
bin_id <- 1:nrow(meth)
strand <- rep(".", nrow(meth))
grange_cpg  <- data.frame(meth[, 1:3], bin_id, CpG, strand)

colnames(grange_cpg)[1] <- "#chr"
write.table(grange_cpg, file = paste0("meth_qc_quant_spikein/", sample_id, "_Granges_CpGs.bed"),
            row.names = F, col.names = T, sep = "\t", quote = F)

## export each quantification column as a separate file
for(i in 4:5) {
    out <- meth[, i]
    feature <- colnames(meth)[i]
    write.table(out, file = paste0("meth_qc_quant_spikein/", sample_id, "_", feature, ".txt"),
                row.names = F, col.names = sample_id, quote = F)
}

# print(object.size(meth_quant),units="auto")
save(meth, file = paste0("meth_qc_quant_spikein/", sample_id, "_meth_quant.RData"))

}
