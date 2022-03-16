##########################################################################
##'@title Function takes a MEDIPSset and finds means of bin counts
##########################################################################
#'@description Function takes MEDIPS SET object, re-organize bins into groups of same CpG counts and finds mean bin counts for each group. This is an internal function used by function "MeDEStrand.plotCalibrationCurve()" and function "MeDEStrand.binMethyl()"
##'@param	MSet \linkS4class{MEDIPSset} objects created by function 'MeDEStrand.createSet()'.
##'@param CSet A \linkS4class{COUPLINGset} object.
##'@return means of bin counts for each group (categorized by bin CpG counts); estimated upper asymptote for the sigmoidal logistic regression model that estimates CpG bias.
##Modified:	July 2017

MeDEStrand.calibrationCurve <- function(MSet=NULL, CSet=NULL, input=FALSE){

	signal=		genome_count(MSet)
	coupling=	genome_CF(CSet)
    maxCoup = floor(max(coupling)*1)  #

	cat("Calculating calibration curve...\n")
	mean_signal = NULL
	coupling_level = NULL

	count_decrease_steps = 0
	max_signal_index = NULL

	first = TRUE
	n_coupling = 0

	for(i in 1:(maxCoup+1)){

		#Test if coupling level is non-empty (can happen for ROIs)
		if(length(signal[coupling==i-1])>0){

			n_coupling = n_coupling + 1
            mean_signal = c(mean_signal, mean(signal[coupling==i-1])) #
			coupling_level = c(coupling_level, i-1)

			##Test if mean_signal decreases
			if(!first){

				if(is.null(max_signal_index) & mean_signal[n_coupling]<mean_signal[n_coupling-1]){
					count_decrease_steps=count_decrease_steps+1
					if(count_decrease_steps==3){
						max_signal_index=n_coupling-3
					}
				}
				else if(is.null(max_signal_index) & mean_signal[n_coupling]>=mean_signal[n_coupling-1]){
					count_decrease_steps=0
				}
			}
			else{
				first = FALSE
			}
		}
# 		else{
#             # cat(paste("Skipping coupling level ", i, " (no data).\n", sep=""))
# 		}
	}

	if(!input & is.null(max_signal_index)){stop("The dependency of coverage signals on sequence pattern (e.g. CpG) densities is different than expected. No Sigmoid model can be build, please check the calibration plot.")}

 		gc()

	return(list(mean_signal=mean_signal, coupling_level=coupling_level, max_index=max_signal_index))
}
