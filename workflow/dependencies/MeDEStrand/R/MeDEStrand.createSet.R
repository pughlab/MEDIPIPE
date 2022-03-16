##########################################################################
##'@title Function to create a list of two MEDIPS SET objects for reads mapped to the positive and negative DNA strand
##########################################################################
##'@description Function creates \linkS4class{MEDIPSset} objects from input data (i.e. aligned short reads) mapped to the positive and negative DNA strand respectively.
##'@description Input format:	bam file or tab (|) separated txt file "chr | start | stop  | strand"
##'@param	file Path and file name of the input data.
##'@param extend defines the number of bases by which the region will be extended before the genome vector is calculated.Regions will be extended along the plus or the minus strand as defined by their provided strand information.
##'@param shift As an alternative to the extend parameter, the shift parameter can be specified. Here, the reads are not extended but shifted by the specified number of nucleotides with respect to the given strand infomation.One of the two parameters extend or shift has to be 0.
##'@param window_size i.e. bin size. This parameter defines the genomic resolution by which short read coverage is calculated.
##'@param BSgenome The reference genome name as defined by BSgenome.
##'@param uniq The uniq parameter determines, if all reads mapping to exactly the same genomic position should be kept (uniq = 0), replaced by only one representative (uniq = 1), or if the number of stacked reads should be capped by a maximal number of stacked reads per genomic position determined by a poisson distribution of stacked reads genome wide and by a given p-value (1 > uniq > 0) (deafult: 1e-3). The smaller the p-value, the more reads at the same genomic position are potentially allowed.
##'@param chr.select only data at the specified chromosomes will be processed.
##'@param paired option for paired end reads.
##'@param sample_name name of the sample to be stored with the MEDIPS SET.
##'@param isSecondaryAlignment option to import only primary alignments.
##'@param simpleCigar option to import only alignments with simple Cigar string.
##'@return \linkS4class{MEDIPSset} objects created for reads mapped to the positive and negative DNA strand respectively
##'@examples # Set parameters
#'@examples file.name = "ENCFF002BKV.bam"; BSgenome="BSgenome.Hsapiens.UCSC.hg19"; uniq=1; extend=200; shift=0; ws=100;
#'chr.select="chr10"
##'@examples # Generate a list of two MEDIP SET objects from reads mapped to the positive and negative DNA strand
##',respectively.
#'@examples MeDIPSset = MeDEStrand.createSet(file=file.name, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq,
#'window_size=ws, chr.select=chr.select, paired = F)
## Requires:	gtools, BSgenome
## Modified:	July 2017


MeDEStrand.createSet <-
function (file = NULL, extend = 0, shift = 0, window_size = 300,
    BSgenome = NULL, uniq = 1e-3, chr.select = NULL, paired = F,
    sample_name = NULL, isSecondaryAlignment = FALSE, simpleCigar=TRUE)
{
    if (is.null(BSgenome)) {
        stop("Must specify a BSgenome library.")
    }
    if (is.null(file)) {
        stop("Must specify a bam or txt file.")
    }
    fileName = unlist(strsplit(file, "/"))[length(unlist(strsplit(file,
        "/")))]
    path = paste(unlist(strsplit(file, "/"))[1:(length(unlist(strsplit(file,
        "/")))) - 1], collapse = "/")
    if (path == "") {
        path = getwd()
    }
    ## dataset = get(ls(paste("package:", BSgenome, sep = "")))
    ## changed to
    dataset = get(ls(paste("package:", BSgenome, sep = ""))[1])
    if (is.null(chr.select)) {
        cat("All chromosomes in the reference BSgenome will be processed:\n")
        chr.select = seqnames(dataset)
        print(chr.select)
    }
    else {
        if (sum(!chr.select %in% seqnames(dataset)) != 0) {
            cat("The requested chromosome(s)", chr.select[!chr.select %in%
                seqnames(dataset)], "are not available in the BSgenome reference. Chromosomes available in the BSgenome reference:",
                seqnames(dataset), "\n")
            stop("Please check chromosome names.")
        }
    }
    if (length(chr.select) > 1) {
        chr.select = gtools::mixedsort(unique(chr.select))
    }
    if (!fileName %in% dir(path)) {
        stop(paste("File", fileName, " not found in", path, sep = " "))
    }
    ext = strsplit(fileName, ".", fixed = T)[[1]]
    if (ext[length(ext)] %in% c("gz", "zip", "bz2"))
        ext = ext[-length(ext)]
    if (ext[length(ext)] %in% c("wig", "bw", "bigwig")) {
        MEDIPSsetObj = getMObjectFromWIG(fileName, path, chr.select,
            BSgenome)
    }
    else {
        if (!paired) {
            GRange.Reads = MeDEStrand.getGRange(fileName, path, extend,
                shift, chr.select, dataset, uniq, isSecondaryAlignment = isSecondaryAlignment, simpleCigar=simpleCigar)

                pos.GRange.Reads = GRange.Reads[strand(GRange.Reads) =='+', ] # modified

                neg.GRange.Reads = GRange.Reads[strand(GRange.Reads) =='-', ] # modified


        }
        else {
            GRange.Reads = MeDEStrand.getPairedGRange(fileName, path, extend,
                shift, chr.select, dataset, uniq, isSecondaryAlignment = isSecondaryAlignment, simpleCigar=simpleCigar)

                pos.GRange.Reads = GRange.Reads[strand(GRange.Reads) =='+', ] # modified

                neg.GRange.Reads = GRange.Reads[strand(GRange.Reads) =='-', ] # modified


        }
        chr_lengths = as.numeric(seqlengths(dataset)[chr.select])
        no_chr_windows = ceiling(chr_lengths/window_size)
        supersize_chr = cumsum(no_chr_windows)
        Granges.genomeVec = MEDIPS.GenomicCoordinates(supersize_chr,
            no_chr_windows, chr.select, chr_lengths, window_size)


       ######

        cat("Calculating short read coverage at genome wide windows for positive strand...\n")
        pos.overlap = countOverlaps(Granges.genomeVec, pos.GRange.Reads)
        pos.datachr = unique(seqlevels(pos.GRange.Reads))

        if (sum(!chr.select %in% pos.datachr) != 0) {
            cat("Please note, no data in the alignment file for chromosome(s):",
                chr.select[!chr.select %in% pos.datachr], "\n")
        }

        if (is.null(sample_name)) {
            sample_name = fileName
        }

        pos.MEDIPSsetObj = new("MEDIPSset", sample_name = sample_name,
            path_name = path, genome_name = BSgenome, number_regions = length(pos.GRange.Reads),
            chr_names = chr.select, chr_lengths = chr_lengths,
            genome_count = pos.overlap, extend = extend, shifted = shift,
            window_size = window_size, uniq = uniq)

       ######

       cat("Calculating short read coverage at genome wide windows for negative strand...\n")
       neg.overlap = countOverlaps(Granges.genomeVec, neg.GRange.Reads)
       neg.datachr = unique(seqlevels(neg.GRange.Reads))

       if (sum(!chr.select %in% neg.datachr) != 0) {
           cat("Please note, no data in the alignment file for chromosome(s):",
           chr.select[!chr.select %in% neg.datachr], "\n")
       }

       if (is.null(sample_name)) {
           sample_name = fileName
       }

       neg.MEDIPSsetObj = new("MEDIPSset", sample_name = sample_name,
       path_name = path, genome_name = BSgenome, number_regions = length(neg.GRange.Reads),
       chr_names = chr.select, chr_lengths = chr_lengths,
       genome_count = neg.overlap, extend = extend, shifted = shift,
       window_size = window_size, uniq = uniq)

       ######


    }

    # return(MEDIPSsetObj)

    return(  c( pos.MEDIPSsetObj, neg.MEDIPSsetObj )    )


}
