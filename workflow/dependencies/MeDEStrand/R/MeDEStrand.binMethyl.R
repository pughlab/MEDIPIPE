##########################################################################
#'@title Function to infer bin-based absolute methylation levels from enrichment signals
##########################################################################
#'@description Function normalizes MeDIP-seq signals by estimating and subsequently correcting CpG bias from the raw signals (i.e. observed bin counts.)
##'@param MSetInput \linkS4class{MEDIPSset} objects created by function 'MeDEStrand.createSet()', reads mapped to the positive and negative DNA strand are processed separately.
##'@param CSet A \linkS4class{COUPLINGset} object.
##'@param ccObj Default NULL. Return of internal function call to \code{\link{MeDEStrand.calibrationCurve}}
##'@param Granges Default FALSE. Return a vector of inferred absolute methylation levels at user-specified bin size. If TRUE, return a \linkS4class{GRanges} object with each bin's genomic coordinate and its absolute methylation levels.
#'@return Absolute methylation levels in a vector
#'@examples file.name = "ENCFF002BKV.bam"; BSgenome="BSgenome.Hsapiens.UCSC.hg19"; uniq=1; extend=200; shift=0; ws=100;
#'chr.select="chr10"
#'@examples MeDIPSset = MeDEStrand.createSet(file=file.name, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq,
#'window_size=ws, chr.select=chr.select, paired = F)
#'@examples Bin_methylation = MeDEStrand.binMethyl(MSet = MeDIPSset, CSet = CS )
## Modified:	July 2017

MeDEStrand.binMethyl <- function(MSetInput=NULL, CSet=NULL, ccObj=NULL, Granges = FALSE){

  for ( i in 1:2 ) {


    if(is.list(MSetInput)){

        MSet=MSetInput[[i]]

    }

	signal =  genome_count(MSet)
	coupling = genome_CF(CSet)


		##Calculate calibration curve
		#####################################
		ccObj = MeDEStrand.calibrationCurve(MSet=MSet, CSet=CSet, input=F)

     ####################   Fitting a simple logistic regression to get initial parameters:

    index.max = which(ccObj$mean_signal== max( ccObj$mean_signal[1:ccObj$max_index] ) )

     MS = ccObj$mean_signal[1:index.max]

     CF = ccObj$coupling_level[1:index.max]

     model.data = data.frame( model.MS =  MS/max( MS), model.CF = CF  )

     logistic.fit = glm(  model.MS ~ model.CF , family=binomial(logit) , data = model.data)




	##Weight signals by linear regression obtained parameters
	####################
     if ( i == 1) { cat("Estimating and correcting CG bias for reads mapped to the DNA positive strand...\n") }

     if ( i == 2) { cat("Estimating and correcting CG bias for reads mapped to the DNA negative strand...\n") }

        estim=numeric(length(ccObj$mean_signal))  # all 0's

	##For the low range of the calibration curve (< index.max) divide counts by observed mean count
	low_range=1:index.max

    estim[low_range]=ccObj$mean_signal[low_range]

	##For the higher range of the calibration curve (>= index.max) divide counts by estimated count
    high_range = ( length(low_range)+1 ):length(estim)

    y.predict = predict( logistic.fit, data.frame( model.CF = ccObj$coupling_level[high_range] ), type ="response"  )*ccObj$mean_signal[ccObj$max_index]

    estim[high_range] = y.predict

    ###################################################

	#rms normalization:
	signal=signal/estim[coupling+1]
	signal[coupling==0]=0


	#Transform the resulting data range into the consistent interval [0:1] by trimming the upper 0.1 quantile of the data
	#######################
	signal = log2(signal)

	signal[is.na(signal)] = 0

	##Shift signals into positive value range
	####################
	minsignal=min(signal[signal!=-Inf])
	signal[signal!=-Inf]=signal[signal!=-Inf]+abs(minsignal)

	##Transform values into the interval [0:1] by compressing the top 0.05% of the signals
	####################
	maxsignal = quantile(signal[signal!=Inf], 0.9995  )
	signal[signal!=Inf & signal>maxsignal]=maxsignal
	signal=round((signal/maxsignal), digits=2)

	##Eliminate -Inf & Inf --> 0
	#######################
	signal[signal==-Inf | signal ==Inf]=0

    if ( i == 1) {  pos.signal = signal }

    if ( i == 2) { neg.signal = signal }

      }


   merged.signal = (pos.signal+neg.signal)/2

   if( !Granges ) {
   return(  merged.signal)}else{

     chr.select = MSet@chr_names
     window_size = window_size(MSet)
     chr_lengths=unname( seqlengths(BSgenome.Hsapiens.UCSC.hg19)[ seqnames(BSgenome.Hsapiens.UCSC.hg19@seqinfo)%in%chr.select ] )
     no_chr_windows = ceiling(chr_lengths/window_size)
     supersize_chr = cumsum(no_chr_windows)
     chromosomes=chr.select
     all.Granges.genomeVec = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, window_size)


     all.Granges.genomeVec$CF = CS@genome_CF

     all.Granges.genomeVec$binMethyl= merged.signal


     return( all.Granges.genomeVec )


   }



}
