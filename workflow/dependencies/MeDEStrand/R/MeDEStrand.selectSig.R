##########################################################################
##'@title Selects putative DMRs (Differential Methylated Regions) from the result table returned by function 'MeDEStrand.meth()'
##########################################################################
##'@description Based on the results table returned by function 'MeDEStrand.meth()', this function selects bins which show significant differential methylation between two groups of MEDIPS SETs. Selection of significant bins follows according to the specified parameters.
##'@param results results table returned by function \code{\link{MeDEStrand.meth}}
##'@param p.value threshold for significant differential methylation. default 0.1
##'@param adj TRUE or FALSE. Whether multiple test correction is conducted based on the methods: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none.
##'@param bg.counts bin filtering parameter. The parameter requires a minimal number of reads for each of the MEDIPS SET groups. To apply this condition, mean of the bin counts per group is considered.
##'@param merge.within.distance  merges significant differential methylated bins within certain bp distance. default is NULL (do not merge).
##'@return result table of bins that are significantly differentially methylated between two groups MSet1 and MSet2.
##'@examples result = MeDEStrand.meth(MSet1=NSCLC_N.MeDIP,  MSet2=NSCLC_T.MeDIP, CSet=CS, chr=c('chr20'),
##'p.adj="fdr", minRowSum= 12)
##'@examples result.sig = MeDEStrand.selectSig(results = results, p.value = 0.1, adj = T, ratio = NULL, bg.counts = 1,
##'merge.within.distance = NULL )
##Modified:	July 2017


MeDEStrand.selectSig = function(results=NULL, p.value=0.1, adj=T, ratio=NULL, bg.counts=NULL, merge.within.distance = NULL ){

	if(is.null(results)){stop("Must specify a result table.")}
	cat(paste("Total number of windows: ", nrow(results), "\n", sep=""))

	##Background counts- preparation
        ################################
        if(!is.null(bg.counts)){
                if(is.numeric(bg.counts)){
                        bg.t = bg.counts
                }

        }

	## Filter for windows tested for differential  methylation
	##########################################################

	results=results[!is.na(results[,grep("p.value", colnames(results))[1]]),]

	cat(paste("Number of windows tested for differential methylation: ", nrow(results), "\n", sep=""))


	# Filter for p.values
	#######################

	if(adj){value_pvalue = grep("adj.p.value", colnames(results))}else{value_pvalue = grep("p.value", colnames(results))[1]}

	results = results[results[,value_pvalue]<= p.value,]

	if(adj){cat(paste("Remaining number of windows with adjusted p.value<=", p.value, ": ", nrow(results), "\n", sep=""))}else{cat(paste("Remaining number of windows with p.value<=", p.value, ": ", nrow(results), "\n", sep=""))}

	## Ratio  filter
	##############################
	## Filter for ratio
	if( !is.null(ratio) ){

		column_ratio = grep("score.log2.ratio", colnames(results))
		if(length(column_ratio)==0){
			column_ratio = grep("edgeR.logFC", colnames(results))
		}

		results=results[results[,column_ratio]>=log2(ratio) | results[,column_ratio]<=(log2(1/ratio)),]
		cat(paste("Remaining number of windows with ratio >=",ratio," (or <=", round(1/ratio, digits=2), ", respectively): ", nrow(results), "\n", sep=""))
	}


	##  Background counts
  ########################################
        if(!is.null(bg.counts)){
                results=results[results$MSets1.counts.mean>=bg.t | results$MSets2.counts.mean>=bg.t,]
                cat(paste("Remaining number of windows where the mean count of at least one group is >=", round(bg.t, digits=2), ": ", nrow(results),"\n", sep=""))

                }

	gc()

	#####


	if(!is.null(merge.within.distance)){

	  if( nrow(results) ==0 ) { stop( "Bin merge cannot be performed: No differential bins are found. Please loose the filter conditions ( use larger p.value or/and smaller ratio or/and smaller bg.counts )"    )   }else{

	     results = MEDIPS.mergeFrames(frames = results, distance = merge.within.distance  )

	  }
	}

	return(as.data.frame(results,  stringsAsFactors=F))

}


