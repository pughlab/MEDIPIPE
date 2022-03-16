#######################################
##Read bed file
#######################################
##Input:	tab (|) separated bed file "chr | start | stop | name | score | strand | ..."
##Param:	allignment.file, path, extend, shift, uniq, dataset
##Output:	Granges object
##Requires:	GenomicRanges
##Modified:	July 2017


MeDEStrand.getGRange <-
function (fileName, path = NULL, extend, shift, chr.select = NULL, 
    dataset = NULL, uniq = 1e-3, ROI = NULL, isSecondaryAlignment = FALSE, simpleCigar=TRUE) 
{
    ext = substr(fileName, nchar(fileName) - 3, nchar(fileName))
    bam = (ext == ".bam" | ext == ".BAM")
    bamindex = bam & file.exists(paste(path, "/", fileName, ".bai", 
        sep = ""))
    if (bam) {
        scanFlag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = isSecondaryAlignment)
        if (bamindex & (!is.null(chr.select) | !is.null(ROI))) {
            if (!is.null(ROI)) {
                cat("Reading bam alignment", fileName, "\n considering ROIs using bam index\n")
                if (!is.null(extend)) {
                  ROI[, 2] = ROI[, 2] - extend
                  ROI[, 3] = ROI[, 3] + extend
                }
                if (!is.null(shift)) {
                  ROI[, 2] = ROI[, 2] - shift
                  ROI[, 3] = ROI[, 3] - shift
                }
                sel = GRanges(chr.select, IRanges(1, 536870912))
            }
            else {
                cat("Reading bam alignment", fileName, "\n considering ", 
                  chr.select, " using bam index\n")
                sel = GRanges(chr.select, IRanges(1, 536870912))
            }
            scanParam = ScanBamParam(flag = scanFlag, simpleCigar= simpleCigar, what = c("rname", 
                "pos", "strand", "qwidth", "isize", "mpos"), which = sel)
        }
        else {
            cat("Reading bam alignment", fileName, "\n")
            scanParam = ScanBamParam(flag = scanFlag, simpleCigar= simpleCigar, what = c("rname", 
                "pos", "strand", "qwidth", "isize", "mpos"))
        }
        regions = scanBam(file = paste(path, fileName, sep = "/"), 
            param = scanParam)
        regions = do.call(rbind, lapply(regions, as.data.frame, 
            stringsAsFactors = F))
        regions = data.frame(chr = as.character(as.vector(regions$rname)), 
            start = as.numeric(as.vector(regions$pos)), stop = as.numeric(as.vector(regions$pos) + 
                as.vector(regions$qwidth) - 1), strand = as.character(as.vector(regions$strand)), 
            stringsAsFactors = F)
    }
    else {
        cat("Reading bed alignment", fileName, "\n")
        regions = read.table(paste(path, fileName, sep = "/"), 
            sep = "\t", header = FALSE, row.names = NULL, comment.char = "", 
            colClasses = c("character", "numeric", "numeric", 
                "NULL", "NULL", "character"))
        names(regions) = c("chr", "start", "stop", "strand")
    }
    if (!is.null(chr.select) & !bamindex) {
        cat("Selecting ", chr.select, "\n")
        regions = regions[regions[, 1] %in% as.vector(chr.select), 
            ]
    }
    cat("Total number of imported short reads: ", nrow(regions), 
        "\n", sep = "")
    regions = adjustReads(regions, extend, shift)
    cat("Creating GRange Object...\n")
    regions_GRange = GRanges(seqnames = regions$chr, ranges = IRanges(start = regions$start, 
        end = regions$stop), strand = regions$strand)
   
    if(is.logical(uniq)){stop("Parameter 'uniq' is not logical anymore, please specify a p-value and see the MeDEStrand vignette.")}
    if (uniq == 1) {
		cat("Keep at most one 1 read mapping to the same genomic location.\n", sep = "")
		regions_GRange = unique(regions_GRange)
		cat("Number of remaining reads: ", length(regions_GRange), 
			"\n", sep = "")
	} else if (uniq < 1 & uniq > 0) {
		max_dup_number = qpois(1 - as.numeric(uniq), length(regions_GRange) / 
			sum(as.numeric(seqlengths(dataset)[chr.select])))
		max_dup_number = max(1, max_dup_number)
		cat("Keep at most ", max_dup_number, 
			" read(s) mapping to the same genomic location\n", sep = "")
		uniq_regions = unique(regions_GRange)
		dup_number = countMatches(uniq_regions, regions_GRange)
		dup_number[dup_number > max_dup_number] = max_dup_number
		regions_GRange = rep(uniq_regions, times = dup_number)
		cat("Number of remaining reads: ", length(regions_GRange), 
			"\n", sep = "")
	} else if (uniq == 0) {
		cat("Do not correct for potential PCR artefacts (keep all reads).\n", sep = "")
	} else {
		stop("Must specify a valid value for parameter uniq. Please check MeDEStrand vignette.")
	}
    #strand(regions_GRange) = "*"
	return(regions_GRange)
}

MeDEStrand.getPairedGRange <-
function (fileName, path = NULL, extend, shift, chr.select = NULL, 
    dataset = NULL, uniq = 1e-3, ROI = NULL, isSecondaryAlignment = FALSE, simpleCigar=TRUE) 
{
    ext = substr(fileName, nchar(fileName) - 3, nchar(fileName))
    bam = (ext == ".bam" | ext == ".BAM")
    bamindex = bam & file.exists(paste(path, "/", fileName, ".bai", 
        sep = ""))
    if (bam) {
        scanFlag = scanBamFlag(isPaired = T, isProperPair = TRUE, 
            hasUnmappedMate = FALSE, isUnmappedQuery = F, isFirstMateRead = T, 
            isSecondMateRead = F, isSecondaryAlignment = isSecondaryAlignment)
        if (bamindex & (!is.null(chr.select) | !is.null(ROI))) {
            if (!is.null(ROI)) {
                cat("Reading bam alignment", fileName, "\n considering ROIs using bam index\n")
                if (!is.null(extend)) {
                  ROI[, 2] = ROI[, 2] - extend
                  ROI[, 3] = ROI[, 3] + extend
                }
                if (!is.null(shift)) {
                  ROI[, 2] = ROI[, 2] - shift
                  ROI[, 3] = ROI[, 3] - shift
                }
                sel = GRanges(ROI[, 1], IRanges(start = ROI[, 
                  2], end = ROI[, 3]))
            }
            else {
                cat("Reading bam alignment", fileName, "\n considering ", 
                  chr.select, " using bam index\n")
                sel = GRanges(chr.select, IRanges(1, 536870912))
            }
        scanParam = ScanBamParam(flag = scanFlag, simpleCigar= simpleCigar, what = c("rname", 
		    "pos", "strand", "qwidth", "isize", "mpos"), which = sel)    
        }
        else {
            cat("Reading bam alignment", fileName, "\n")
            scanParam = ScanBamParam(flag = scanFlag, simpleCigar = simpleCigar, what = c("rname", 
                  "pos", "strand", "qwidth", "isize", "mpos"))
        }
   
        regions = scanBam(file = paste(path, fileName, sep = "/"), 
            param = scanParam)
        regions = do.call(rbind, lapply(regions, as.data.frame, 
            stringsAsFactors = F))
    }
    else {
        stop("BED files in paired end mode not supported.\n")
    }
    if (!is.null(chr.select) & !bamindex) {
        cat("Selecting", chr.select, "\n")
        regions = regions[regions[, 1] %in% as.vector(chr.select), 
            ]
    }
    cat("Total number of imported first mate reads in properly mapped pairs: ", 
        nrow(regions), "\n", sep = "")
    cat("scanBamFlag: isPaired = T, isProperPair=TRUE , hasUnmappedMate=FALSE, ", 
        "isUnmappedQuery = F, isFirstMateRead = T, isSecondMateRead = F\n", 
        sep = "")
    cat("Mean insertion size: ", mean(abs(regions$isize)), " nt\n", 
        sep = "")
    cat("SD of the insertion size: ", sd(abs(regions$isize)), 
        " nt\n", sep = "")
    cat("Max insertion size: ", max(abs(regions$isize)), " nt\n", 
        sep = "")
    cat("Min insertion size: ", min(abs(regions$isize)), " nt\n", 
        sep = "")
   
   qwidth = regions[, "qwidth"]
   regions = data.frame(chr = as.character(as.vector(regions$rname)), 
   start = as.numeric(as.vector(regions$pos)), stop = as.numeric(as.vector(regions$mpos)), 
   strand = as.character(as.vector(regions$strand)), 
   isize = as.numeric(as.vector(regions$isize)), stringsAsFactors = F)
   
   regionsToRev = regions$start > regions$stop
   regions[regionsToRev, ]$start = regions[regionsToRev,]$stop
   regions[, "stop"] = regions[, "start"] + abs(regions[, "isize"]) - 1    
   
   if(extend!=0){cat("The extend parameter will be neglected, because the actual DNA fragment length is known in paired-end data.\n")}
   if(shift!=0){cat("The shift parameter will be neglected, because the actual DNA fragment position is known in paired-end data.\n")}
    
    cat("Creating GRange Object...\n")
    regions_GRange = GRanges(seqnames = regions$chr, ranges = IRanges(start = regions$start, 
        end = regions$stop), strand = regions$strand)
	
	if(is.logical(uniq)){stop("Parameter 'uniq' is not logical anymore, please specify a p-value and see the MeDEStrand vignette.")}
	
	if (uniq == 1) {
		cat("Keep at most 1 read mapping to the same genomic location.\n", sep = "")
		regions_GRange = unique(regions_GRange)
		cat("Number of remaining short reads: ", length(regions_GRange), 
			"\n", sep = "")
	} else if (uniq < 1 & uniq > 0) {
		max_dup_number = qpois(1 - as.numeric(uniq), length(regions_GRange) / 
			sum(as.numeric(seqlengths(dataset)[chr.select])))
		max_dup_number = max(1, max_dup_number)
		cat("Keep at most ", max_dup_number, 
			" first mate read(s) mapping to the same genomic location\n", sep = "")		
		uniq_regions = unique(regions_GRange)
		dup_number = countMatches(uniq_regions, regions_GRange)
		dup_number[dup_number > max_dup_number] = max_dup_number
		regions_GRange = rep(uniq_regions, times = dup_number)
		cat("Number of remaining short reads: ", length(regions_GRange), 
			"\n", sep = "")
	} else if (uniq == 0) {
		cat("Do not correct for potential PCR artefacts (keep all reads).\n", sep = "")
	} else {
		stop("Must specify a valid value for parameter uniq. Please check MeDEStrand vignette.")
	}
    #strand(regions_GRange) = "*"
	return(regions_GRange)
}

scanBamToGRanges <- function(...) {
	dat <- scanBam(...)[[1]]
	keep <- !is.na(dat$pos)
	GRanges(seqnames=dat$rname[keep],
	ranges=IRanges(start=dat$pos[keep], width=nchar(dat$seq[keep])),
	strand=dat$strand[keep], isize=dat$isize[keep],
	mrnm=dat$mrnm[keep], flag=dat$flag[keep])
}

