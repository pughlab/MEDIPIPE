##########################################################################
##'@title Function plots the 'calibration plot', which reveals the CpG density bias
##########################################################################
#'@description Visualizes the dependency between MeDIP-seq signals and CpG densities by 'calibration curve'. 'calibration curve' estimates CpG bias by fitting a sigmoidal logistic regression curve.
##'@param MSet \linkS4class{MEDIPSset} objects created by function 'MeDEStrand.createSet()' for reads mapped to the positive and negatie DNA strand respectively.
##'@param CSet A \linkS4class{COUPLINGset} object. If NULL, it can be calculated and supplied from MSet.
##'@param Forward_strand if TRUE, plot the 'calibration plot' for the reads mapped to the positive DNA strand(i.e. estimating CpG density bias for reads mapped to the positive DNA strand); if FALSE, plot the 'calibration plot' for the reads mapped to the negative DNA strand (i.e. estimating CpG density bias for reads mapped to the negative DNA strand).
##'@param default = "all". It is recommended to call a graphics device (e.g. png("calibrationPlot.png")) before calling the plot command, because R might not be able to plot the full amount of data in reasonable time.
##'@param main The title of the calibration plot.
##'@param xrange The signal range of the calibration curve typically falls into a low signal range. By setting the xrange parameter to TRUE, the calibration plot will visualize the low signal range only.
##'@return Calibration plot, i.e. Visualization of the CpG bias curve.
#'@examples file.name = "ENCFF002BKV.bam"; BSgenome="BSgenome.Hsapiens.UCSC.hg19"; uniq=1; extend=200; shift=0; ws=100;
#'chr.select="chr10"
#'@examples MeDIPSset = MeDEStrand.createSet(file=file.name, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq,
#'window_size=ws, chr.select=chr.select, paired = F)
#'@examples # plot the 'calibration plot'
#'@examples MeDEStrand.plotCalibrationCurve( MSet=MeDIPset, CSet=NULL,  Forward_Strand = T, main = NULL, xrange=TRUE)
##Modified:	July 2017

MeDEStrand.plotCalibrationCurve <- function(MSet=NULL, CSet=NULL, Forward_Strand = T, plot_chr="all", main="Calibration plot", xrange=T){

  ##Proof data accordance....
  ##########################
  input=F


    if( is.list( MSet ) & ( Forward_Strand == T) ){

      MSet=MSet[[1]]

    }else{ MSet=MSet[[2]]   }

  if(!(is.null(MSet)) & (class(MSet)!="MEDIPSset" )) stop("MSet must be a MEDIPSset object!")

  if(is.null(MSet) ) stop("MSet must be specified")

  #print("Checking CSet")
  if(class(CSet)!="COUPLINGset") stop("Must specify a COUPLINGset object!")

  #print("Checking MSet")
  if(!is.null(MSet) & class(MSet)=="MEDIPSset"){
    if(window_size(MSet)!=window_size(CSet)) stop("MSet and COUPLINGset have different window sizes!")
    if(length(chr_names(MSet))!=length(chr_names(CSet))) stop("MSet and COUPLINGset contain a different number of chromosomes!")
    for(i in 1:length(chr_names(MSet))){
      if(chr_names(MSet)[i]!=chr_names(CSet)[i]){stop("MSset and COUPLINGset contain different chomosomes!")}
    }
  }



  if(!is.null(MSet)){
    signal=genome_count(MSet)
    chr_lengths = chr_lengths(MSet)
    if(class(MSet)=="MEDIPSset"){window_size = window_size(MSet)}
    number_regions = number_regions(MSet)
    chromosomes=chr_names(MSet)
  }
  coupling=genome_CF(CSet)
  seq_pattern=seq_pattern(CSet)

  ##Calculate calibration curve
  #####################################
  if (!is.null(MSet))
    ccObj_MSet = MEDIPS.calibrationCurve(MSet=MSet, CSet=CSet, input=F)


  ##Check, if a subset of chromosomes has been selected
  ######################################################
  if(plot_chr!="all" & (class(MSet)=="MEDIPSset" )){
    cat(paste("Extracting data for",plot_chr, "...\n", sep=" "))

    ##Calculate genomic coordinates
    ##################################
    no_chr_windows = ceiling(chr_lengths/window_size)
    supersize_chr = cumsum(no_chr_windows)
    genome_chr = as.vector(seqnames(MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, window_size)))

    if(length(genome_chr[genome_chr==plot_chr])==0){stop("Stated calibration chromosome does not exist within the MEDIPS SET.")}
    signal = signal[genome_chr==plot_chr]
    coupling = coupling[genome_chr==plot_chr]
    cat(paste("Plotting calibration plot for", plot_chr, "...\n", sep=" "))
  }

  if(plot_chr=="all"){cat("Plotting calibration plot for all chromosomes. It is recommended to redirect the output to a graphic device.\n")}


  descSignal = "# Reads/bin"


  ##Preparations for raw data plots
  #################################
  if(xrange){
    if(!is.null(MSet))
      range=c(0,max(ccObj_MSet$mean_signal)*5)

  }else{
    range=c(0,max(signal))
  }
  ##Plot
  #######
  plot(coupling,signal, pch=".", main=main, ylab=descSignal,ylim=range*0.7, xlab=paste0("bin CpG counts"), col="lightblue" , cex.lab = 1.2)
  for(i in 0:max(coupling)){
    t=table(signal[coupling==i])
    points(x=rep(i,length(t)), y=as.numeric(names(t)),lwd=log(t, 10), pch=4, col="lightblue")
  }
  if(!is.null(MSet))
    llab="MeDIP-seq bin reads"
  else
    llab="Input reads in genomic bin"
  lcol="lightblue"
  if(! is.null(MSet)){

    index.max = which(ccObj_MSet$mean_signal== max( ccObj_MSet$mean_signal[1:ccObj_MSet$max_index] ) )

    MS = ccObj_MSet$mean_signal[1:index.max]

    CF = ccObj_MSet$coupling_level[1:index.max]
    model.data = data.frame( model.MS =  MS/max( MS), model.CF = CF  )

    logistic.fit = glm(  model.MS ~ model.CF , family=binomial(logit) , data = model.data)


    lines(ccObj_MSet$coupling_level, ccObj_MSet$mean_signal, col="red",lwd = 1.5)


    lines( ccObj_MSet$coupling_level , predict( logistic.fit, data.frame( model.CF = ccObj_MSet$coupling_level), type ="response"  )*ccObj_MSet$mean_signal[index.max] ,  col = 'blue' ,lwd =2 )


    llab=c("Mean bin reads","Estimated CpG bias")


    lcol=c("red","blue")


  }

  legend("topright", legend=llab, fill=lcol)

}
