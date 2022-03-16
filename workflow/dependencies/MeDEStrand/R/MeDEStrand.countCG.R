##########################################################################
##'@title The function counts the frequency of the specified sequence pattern exist in each bin.
##########################################################################
##'@description The function calculates the local densities(i.e. bin level) of a defined sequence pattern (e.g. 'CG' for CpGs) and returns a \linkS4class{COUPLINGset} object.
##Input:	Parameters that specify the sequence pattern and the size of the bins w.r.t the reference genome
##'@param	pattern defines the sequence pattern, e.g. 'CG' for CpGs.
##'@param refObj A \linkS4class{MEDIPSset} object or an object returned by function \code{\link{MeDEStrand.createSet}}
##'@return 	the counts of specified pattern e.g.'CG' frequencies in bins
##'@examples file.name = "ENCFF002BKV.bam"; BSgenome="BSgenome.Hsapiens.UCSC.hg19"; uniq=1; extend=200; shift=0; ws=100;
##'chr.select="chr10"
##'@examples MeDIPSset = MeDEStrand.createSet(file=file.name, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq,
##'window_size=ws, chr.select=chr.select,paired = F)
##'@examples # count CG contents in each bin
##'@examples CS = MeDEStrand.countCG(pattern="CG", refObj=MeDIPSset)
##Requires:	gtools, BSgenome
##Modified:	July 2017

MeDEStrand.countCG <- function(pattern="CG", refObj=NULL){
	if(is.list(refObj)){
		refObj=refObj[[1]]
	}
	if(class(refObj)!="MEDIPSset")	{
		stop("You must provide an MEDIPSset as reference object\n")
	}

	chr_lengths=chr_lengths(refObj)
	chromosomes = chr_names(refObj)
	BSgenome = genome_name(refObj)

	if(is.list(refObj)){
		refObj=refObj[[1]]
	}

	## Get the genomic positions of the sequence pattern
	cat("Get genomic sequence pattern positions...\n")
	GRanges.pattern = MEDIPS.getPositions(BSgenome, pattern, chromosomes)

	## Create the genome vector Granges object
	if(class(refObj)=="MEDIPSset"){
		window_size = window_size(refObj)
		no_chr_windows = ceiling(chr_lengths/window_size)
		supersize_chr = cumsum(no_chr_windows)
		Granges.genomeVec = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, window_size)
	}else{
       stop("You must provide an MEDIPSset as reference object\n")
	}

	##Count the number of sequence pattern in each window.
	cat(paste("Counting the number of ",pattern, "'s in each window...\n", sep=""))
	genomeCoup = countOverlaps(Granges.genomeVec, GRanges.pattern)

	COUPLINGsetObj = new('COUPLINGset', seq_pattern=pattern, genome_name=BSgenome, genome_CF=genomeCoup, number_pattern=length(GRanges.pattern), window_size=window_size, chr_names=chromosomes, chr_lengths=chr_lengths)
	gc()

	return(COUPLINGsetObj)
}
