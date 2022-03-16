##########################################################################
##'@title Function calculates genome wide absolute methylation levels for each provided set and uses t test for differential methylated loci between two groups of MEDIPS SETs (if provided)
##########################################################################
##'@description The function summarizes coverage profiles (bin counts) as well as infers absolute methylation levels for given MEDIPS SETs. In case two groups of MEDIPS SETs are provided at MSet1 and MSet2, the function finds differential methylated loci between MSet1 and MSet2 based on the inferred absolute methylation levels. For each bin/loci, t test is conducted and multiple test correction ( genome-wide/chromosome-wide ) can be specified by user.
##'@param MSet1 first group of MEDIPS SETs. Each group contains several \linkS4class{MEDIPSset} objects.
##'@param MSet2 second group of MEDIPS SETs. Each group contains several \linkS4class{MEDIPSset} objects. Differential methylated coverage will be calculated by t test and returned, if MSet1 and MSet2 are not empty.
##'@param CSet A \linkS4class{COUPLINGset} object. (i.e. corresponding bin CpG density for the provided MSet1 or MSet2). This object is return of function \code{\link{MeDEStrand.countCG}}
##'@param chr specifies the chromosome(s) for t test for differential methylated loci between MSet1 and MSet2. Depending on chromosome(s) or all chromosomes are provided, multiple test correction are conducted at chromosome-wide or genome-wide levels. The latter one is more stringent.
##'@param p.adj in order to correct p.values for multiple testing, MeDEStrand uses R's 'p.adjust()' function. Therefore, the following methods are available: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none.
##'@param minRowSum threshold for a minimum sum of counts across all samples bins (default=10). Bins with lower coverage will not be tested for differential methylation coverage.
##'@return 	result table of summary statistics and t test results
##'@examples result = MeDEStrand.meth(MSet1=NSCLC_N.MeDIP,  MSet2=NSCLC_T.MeDIP, CSet=CS, chr=c('chr20'),
##'p.adj="fdr", minRowSum= 12)
##Modified:	July 2017

MeDEStrand.meth = function(
  MSet1 = NULL,
  MSet2 = NULL,
  CSet = NULL,

  chr = NULL,
  p.adj="bonferroni",

  minRowSum=10

)
{
  # MSet1 = NSCLC_N.MeDIP; MSet2 = NSCLC_T.MeDIP; CSet = CS ;    chr=c('chr20','chr21')

  nMSets1 = length(MSet1)/2
  nMSets2 = length(MSet2)/2

  if(is.list(CSet)) CSet=CSet[[1]]
  if(!is.list(MSet1)) MSet1=c(MSet1)
  if(!is.list(MSet2)) MSet2=c(MSet2)

  controlSet = MSet1[[1]]

  ##################################
  ##Calculate genomic coordinates
  ##################################

  if(class(controlSet)=="MEDIPSset"){
    window_size = window_size(controlSet)
    no_chr_windows = ceiling(chr_lengths(controlSet)/window_size(controlSet))
    supersize_chr = cumsum(no_chr_windows)
    GRanges.genome = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chr_names(controlSet), chr_lengths(controlSet), window_size(controlSet))
  }

  ##Create data frame for all genomic windows and MEDIPS SETs

  base = data.frame(chr=as.vector(seqnames(GRanges.genome)), start=start(GRanges.genome), end=end(GRanges.genome), CF=genome_CF(CSet), stringsAsFactors=F)

  if(!is.null(chr)){
    fi=base[,1]%in%chr

    cat("Extracting data for", chr, "...\n", sep=" ")
    if(sum(fi)==0){stop("Stated chromosome does not exist in the supplied SET.")}

  }

  counts.medip = NULL
  rms = NULL



  ##Add  inferred absolute methylation levels
  if(!is.null(MSet1)){
    for(i in 1:nMSets1 ){
      cat(paste("Preprocessing MEDIPS SET ", i, " in MSet1...\n", sep=""))

      counts.medip = cbind(counts.medip, MSet1=genome_count(MSet1[[2*i-1]]) + genome_count(MSet1[[2*i]])     )
    #  ccObj = MEDIPS.calibrationCurve(MSet=MSet1[[i]], CSet=CSet, input=F)
      rms = cbind(rms, MeDEStrand.binMethyl( c(MSet1[[2*i-1]],MSet1[[2*i]]  ), CSet ) )

    }
  }

  if(!is.null(MSet2)){
    for(i in 1:nMSets2 ){

      cat(paste("Preprocessing MEDIPS SET ", i, " in MSet2...\n", sep=""))

      counts.medip = cbind(counts.medip, MSet2=genome_count(MSet2[[2*i-1]]) + genome_count(MSet2[[2*i]])  )
   #   ccObj = MEDIPS.calibrationCurve(MSet=MSet2[[i]], CSet=CSet, input=F)
      rms = cbind(rms, MeDEStrand.binMethyl(c(MSet2[[2*i-1]],MSet2[[2*i]]  ), CSet  ))

    }
  }


  ##Extract data for selected chromosome
  #######################################
  if(!is.null(chr)){

    fi=base[,1]%in%chr

    cat("Extracting data for", chr, "...\n", sep=" ")
    if(sum(fi)==0){stop("Stated chromosome does not exist in the supplied SET.")}

    if(!is.null(counts.medip)){

      counts.medip = counts.medip[fi,]
      rms = rms[fi,]
    }

    base = base[fi,]
    cat(nrow(base), "windows on", chr, "\n",sep=" ")
  }

  ##Set colnames and transform to data.frames
  ############################################
  col.names.count = NULL
  col.names.rms = NULL
  if(nMSets1!=0){
    for(i in 1:nMSets1 ){

      col.names.count = c(col.names.count, paste(sample_name(MSet1[[2*i-1]]), ".counts", sep=""))
      col.names.rms = c(col.names.rms, paste(sample_name(MSet1[[2*i-1]]), ".rms", sep=""))

    }
  }

  if(nMSets2!=0){
    for(i in 1:nMSets2 ){

      col.names.count = c(col.names.count, paste(sample_name(MSet2[[2*i-1]]), ".counts", sep=""))
      col.names.rms = c(col.names.rms, paste(sample_name(MSet2[[2*i-1]]), ".rms", sep=""))

    }
  }

  if(nMSets1!=0 | nMSets2!=0){

    counts.medip = data.frame(counts.medip)
    colnames(counts.medip) = col.names.count

    rms = data.frame(rms)
    colnames(rms) = col.names.rms

  }


  ## If two groups of MEDIPS SETs are given
  ## calculate differential coverage
  ##################################
  if(!is.null(MSet1) & !is.null(MSet2)){

    cat(paste("Differential coverage analysis...\n", sep=" "))

    ## Correct for test selection if necessary

    diff.results.list = MeDEStrand.diffMeth(base=base, values= rms, bin.counts= counts.medip  , nMSets1=nMSets1, nMSets2=nMSets2, p.adj=p.adj, minRowSum=minRowSum)

    cat("Please note, log2 ratios are reported as log2(MSet1/MSet2).\n")
    diff.results = diff.results.list$diff.results
    diff.index = diff.results.list$diff.index

  }



  ##Create results table
  ##################################
  cat(paste("Creating results table...\n", sep=" "))
  if(!is.null(counts.medip)){

    results = data.frame(base, counts.medip, rms, stringsAsFactors=F)


  }



  ##Add mean counts, mean rms columns
  set1idx=1:(nMSets1)
  counts.mean.C=numeric(dim(counts.medip)[1])

  for (i in set1idx){
    counts.mean.C=counts.mean.C+counts.medip[,i]

  }
  counts.mean.C=counts.mean.C/nMSets1

  rms.mean.C = numeric(dim(rms)[1])
  for (i in set1idx){
    rms.mean.C=rms.mean.C+rms[,i]
  }
  rms.mean.C=rms.mean.C/nMSets1
  results = data.frame(results, MSets1.counts.mean=counts.mean.C,  MSets1.rms.mean=rms.mean.C, stringsAsFactors=F)



  if(nMSets2>0){
    set2idx=(nMSets1+1):(nMSets1+nMSets2)
    counts.mean.T=numeric(dim(counts.medip)[1])

    for (i in set2idx){
      counts.mean.T=counts.mean.T+counts.medip[,i]

    }
    counts.mean.T=counts.mean.T/nMSets2

    rms.mean.T =numeric(dim(rms)[1])
    for (i in set2idx){
      rms.mean.T=rms.mean.T+rms[,i]
    }
    rms.mean.T=rms.mean.T/nMSets2
    results = data.frame(results, MSets2.counts.mean=counts.mean.T,  MSets2.rms.mean=rms.mean.T, stringsAsFactors=F)


  }


################ Good up here


  ##Add diff.meth results
  if(nMSets1!=0 & nMSets2!=0){
    cat(paste("Adding differential coverage results...\n", sep=" "))
    dummy.results = matrix(ncol=ncol(diff.results), nrow=nrow(results))


    dummy.results[diff.index,] = diff.results
    colnames(dummy.results)=colnames(diff.results)
    results = data.frame(results, dummy.results, stringsAsFactors=F)

  }



  rownames(results) = seq(1, nrow(results))



  return(results)

}

