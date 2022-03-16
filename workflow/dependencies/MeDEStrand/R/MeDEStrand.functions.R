##########
#Function
##########
##Input:	Parameters that specify the reference genome and targeted window size.
##Param:	supersize_chr, no_chr_windows_, chromosomes, chr_lengths, window_size
##Output:	GRange object
##Modified:	July 2017

matNnotNA <- function(x,dimension=c("row"=1,"col"=2)[1]){
  #returns number of non-na-values per row/col
  return(apply(X=!is.na(x),MARGIN=dimension,FUN=sum))
}

## ---------------------------------------------------
##########
#Function
##########
##Input:	Parameters that specify the reference genome and targeted window size.
##Param:	supersize_chr, no_chr_windows_, chromosomes, chr_lengths, window_size
##Output:	GRange object
##Modified:	July 2017

matMin <- function(x,dimension=c("row"=1,"col"=2)[1]){
  return(apply(X=x,MARGIN=dimension,FUN=min))
}
matDiff<-function(x,groups,dimension=c("row"=1,"col"=2)[1]){
  ugr <- sort(unique(groups))
  f0 <- groups==ugr[1] & !is.na(groups)
  f1 <- groups==ugr[2] & !is.na(groups)
  if(dimension==1){	
    res<-x[,f0,drop=FALSE]-x[,f1,drop=FALSE]
  }else{
    res<-x[f0,,drop=FALSE]-x[f1,,drop=FALSE]
  }
  return(res)
}

## ---------------------------------------------------
matMax <- function(x,dimension=c("row"=1,"col"=2)[1]){
  return(apply(X=x,MARGIN=dimension,FUN=max))
}
## ---------------------------------------------------
##########
#Function
##########
##Input:	Parameters that specify the reference genome and targeted window size.
##Param:	supersize_chr, no_chr_windows_, chromosomes, chr_lengths, window_size
##Output:	GRange object
##Modified:	July 2017

matMean <- function(x,dimension=c("row"=1,"col"=2)[1]){
  return(apply(X=x,MARGIN=dimension,FUN=mean))
}
## ---------------------------------------------------
##########
#Function
##########
##Input:	Parameters that specify the reference genome and targeted window size.
##Param:	supersize_chr, no_chr_windows_, chromosomes, chr_lengths, window_size
##Output:	GRange object
##Modified:	July 2017

matSd <- function(x,dimension=c("row"=1,"col"=2)[1]){
  return(apply(X=x,MARGIN=dimension,FUN=sd))
}
## ---------------------------------------------------
##########
#Function
##########
##Input:	Parameters that specify the reference genome and targeted window size.
##Param:	supersize_chr, no_chr_windows_, chromosomes, chr_lengths, window_size
##Output:	GRange object
##Modified:	July 2017

matTtest <- function(x,groups,dimension=c("row"=1,"col"=2)[1],alternative = c("two.sided", "less", "greater")[1],paired=c(FALSE,TRUE)[1]){
  #cat("modified version!\n")
  if(length(groups)!= dim(x)[c(2,1)[dimension]]){
    stop(paste("Length of groups vector is not identical to dimension",c(2,1)[dimension],"  of x."))
  }
  if(paired){
    if(length(f0)!=length(f1)){
      stop("not all arguments have the same length")
    }
    vdiff<-matDiff(x,groups)
    #######################
    n=apply(!is.na(vdiff),dimension,sum)
    #######################
    tstat=matMean(vdiff)/sqrt(apply(X=x,MARGIN=dimension,FUN=var)/n)
    df <- n -1
  }else{
    ugr <- sort(unique(groups))
    if(length(ugr)!=2){
      stop(paste("There need to be exactly 2 group identifiers in groups (not counting NAs) but",
                 length(ugr),"were provided."))
    }
    f0 <- groups==ugr[1] & !is.na(groups)
    f1 <- groups==ugr[2] & !is.na(groups)
    sd0<- apply(x[,f0],dimension,sd)
    sd1<- apply(x[,f1],dimension,sd)
    na0<-apply(!is.na(x[,f0]),dimension,sum)
    na1=apply(!is.na(x[,f1]),dimension,sum)
    tstat <- (apply(x[,f0],dimension,mean)-apply(x[,f1],dimension,mean)) / sqrt( sd0^2/na0 + sd1^2/na1)
    df <- na0 + na1 -2
    if(any(sd0 != sd1)){
      fi     <- sd0 != sd1 			
      err0   <- sd0[fi]/sqrt(na0[fi])
      err1   <- sd1[fi]/sqrt(na1[fi])		
      df[fi] <- (err0^2+err1^2)^2/(err0^4/(na0[fi]-1) + err1^4/(na1[fi]-1))	
      
    }
  }
  if (alternative == "less") {
    pval <- pt(tstat, df)
  }
  else if (alternative == "greater") {
    pval <- pt(tstat, df, lower.tail = FALSE)
  }
  else {
    pval  <- 2 * pt(-abs(tstat), df)
  }
  return(data.frame(p.value=pval,df=df, t.statistics=tstat))
  ################
  #return(pval)
  
}
