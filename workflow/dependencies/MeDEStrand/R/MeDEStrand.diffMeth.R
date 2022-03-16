##########################################################################
##Function provides two moudiles for calculating differential methylation
##########################################################################
##Input:	genomic coordinates and data for groups of MEDIPS SETs
##Param:	base, values, bin.counts, nMSets1, nMSets2, p.adj, minRowSum
##Output:	results of differential methylation analysis plus index in genome wide table
##Requires:	edgeR
##Modified:	July 2017

MeDEStrand.diffMeth = function(base = NULL, values=NULL, bin.counts= NULL ,  nMSets1=NULL, nMSets2=NULL, p.adj="dfr",  minRowSum=10 )
{


	##ttest/Score##
	#########

    # values = rms; p.adj="fdr";  minRowSum=10

		##Extract non-zero MeDIP rows
		cat(paste("Extracting count windows with at least",minRowSum," reads...\n", sep=" "))

		filter= rowSums(bin.counts)>=minRowSum

		##Extract non-zero coupling factor(CpG) windows

			cat(paste("Extracting non-zero CpG windows...\n", sep=" "))
			filter=filter & base[,4]!=0


		#filter NA
		#if there is a na in a row, this row is sorted out?!
		filter=filter & !is.na(rowSums(values))


		cat(paste("Calculating score for", sum(filter), "windows...\n", sep=" "))
		ms1=1:nMSets1
		ms2=(nMSets1+1):(nMSets1+nMSets2)
		##Calculate ratios##
		####################
		ratio=rowSums(values[filter,ms1]+0.1)/rowSums(values[filter,ms2,drop=F]+0.1)*(nMSets2/nMSets1)

		methylation.change = rep( 1 , length(ratio)  )

		methylation.change[ratio<1 ] = -1


		##Calculate p.values##
		######################
		##Check for constant entries
		## Filter out constant entries, as they have no sd

		const = apply(X=values[filter,ms1,drop=F],MARGIN=1,FUN=min) - apply(X=values[filter,ms1,drop=F],MARGIN=1,FUN=max) 	+
			apply(X=values[filter,ms2,drop=F],MARGIN=1,FUN=min) - apply(X=values[filter,ms2,drop=F],MARGIN=1,FUN=max) 	== 0

		t.test.p.value = rep(NA, sum(filter))

		t.test.p.value[!const]=matTtest(values[filter,][!const,],groups=c(rep(1,nMSets1),rep(2,nMSets2)))$p.value

		##Calculate the final score##
		#############################
		score = (-log10(t.test.p.value)*10)*log(ratio)

		##Adjusting p.values for multiple testing

		cat(paste("Adjusting p.values for multiple testing...\n", sep=" "))

		diff.results = cbind( methylation.change = methylation.change, log2.ratio=log2(ratio), p.value=t.test.p.value, adj.p.value = p.adjust(t.test.p.value, p.adj), score=score)

		rm(const, ratio, t.test.p.value, score)


	return(list(diff.results=diff.results, diff.index=which(filter)))
}

