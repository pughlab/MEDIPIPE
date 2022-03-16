##########################################################################
##'@title Function calculates correlation (Pearson or Spearman) coefficient between inferred bin-based absolute methylation levels with supplied RRBS data
##########################################################################
#'@description Function calculates correlation (Pearson or Spearman) coefficient between a vector of inferred bin-based absolute methylation levels with supplied RRBS data. Since RRBS data provides single-resolution CpG cytosine methylation levels, thus RRBS CpGs are binned and the means of CpG methylation in the bins are assigned as the 'true' absolute methylation levels for the bins.
##'@param MSetInput \linkS4class{MEDIPSset} objects.
##'@param CSet A \linkS4class{COUPLINGset} object.
##'@param RRBS A \linkS4class{methylRaw} object ( object returned by function "processBismarkAln()" from package 'MethylKit') or A \linkS4class{GRanges} object with a metadata column 'methylation' to specify the methylation level of each CpG or a '.bed' file with at least these columns: 'chr','start','end','methylation' for CpGs.
##'@param minRRBS to filter out bins with RRBS CpGs less than (<) 'minRRBS' number.
##'@param method "pearson" or "spearman". To find the Pearson or Spearman correlation coefficient.
#'@return Pearson/Spearman correlation coefficient value
#'@examples # Generate MEDIPS SET objects from .bam file
#'@examples file.name = "ENCFF002BKV.bam"; BSgenome="BSgenome.Hsapiens.UCSC.hg19"; uniq=1; extend=200; shift=0; ws=100;
#'chr.select="chr10"
#'@examples MeDIPSset = MeDEStrand.createSet(file=file.name, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq,
#'window_size=ws, chr.select=chr.select, paired = F)
#'@examples # Load RRBS data (a 'methylRaw' object) saved in the package
#'@examples fpath <- system.file("data", 'GM12878.RRBS.RData', package="MeDEStrand")
#'@examples load(fpath)
#'@examples # Find the Pearson correlation coefficient
#'@examples correlation = MeDEStrand.cor( MSetInput= MeDIP_single, CS = CS, RRBS = RRBS, minRRBS = 4 , method = "pearson" )
## Modified:	July 2017


MeDEStrand.cor <- function( MSetInput=NULL, CS = NULL, RRBS = NULL, minRRBS = NULL , method = "pearson"  )
{

  # Transfer bins to Granged object

  if(is.list(MSetInput)){


    MSet=MSetInput[[1]]

  }


  if( is.null(CS) ) {
  CS = MeDEStrand.countCG(pattern="CG", refObj=MSet) }


  bin.Methylation = MeDEStrand.binMethyl(MSet = MSetInput, CSet = CS )  # original MEDIPS inferred bin methylation


  chr.select = MSet@chr_names


  window_size = window_size(MSet)
  chr_lengths=unname( seqlengths(BSgenome.Hsapiens.UCSC.hg19)[ seqnames(BSgenome.Hsapiens.UCSC.hg19@seqinfo)==chr.select ] )
  no_chr_windows = ceiling(chr_lengths/window_size)
  supersize_chr = cumsum(no_chr_windows)
  chromosomes=chr.select
  all.Granges.genomeVec = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, window_size)
  all.Granges.genomeVec$ID = 1:length(all.Granges.genomeVec)

  correlation = vector()

  for ( CHR in chr.select ) {


  Granges.genomeVec = all.Granges.genomeVec[seqnames(all.Granges.genomeVec) == CHR, ]


  Granges.genomeVec$bin=1:length(Granges.genomeVec)  # the corresponding bin ID

  Granges.genomeVec$CF = CS@genome_CF[Granges.genomeVec$ID]

  Granges.genomeVec$binMethyl=bin.Methylation[Granges.genomeVec$ID]



  ##### Get GRange data type of RRBS data:

  if (class(RRBS)=="methylRaw"){
  G.RRBS = GRanges(seqnames=RRBS$chr, ranges=IRanges(RRBS$start, RRBS$end),
                   strand=RRBS$strand, methylation = RRBS$numCs/RRBS$coverage)

  }

  if (class(RRBS)=="GRanges"){

     if( is.null(RRBS$methylation ))   stop("Please add a metadata column with name 'methylation' for the cytosines' methylation levels betwee 0 and 1")

    }

  if (class(RRBS)=="data.frame"){

    if( is.null(RRBS$chr ))   stop("Please add a metadata column with name 'chr' for the chromosome(s)")
    if( is.null(RRBS$start ))   stop("Please add a metadata column with name 'start' for cytosine position")
    if( is.null(RRBS$end ))   stop("Please add a metadata column with name 'end' for the cytosine position, should be the same as column 'start'")
    if( is.null(RRBS$strand ))   stop("Please add a metadata column with name 'strand' for the cytosines' strand information")
    if( is.null(RRBS$methylation ))   stop("Please add a metadata column with name 'methylation' for the cytosines' methylation levels betwee 0 and 1")

    G.RRBS = GRanges(seqnames=RRBS$chr, ranges=IRanges(RRBS$start, RRBS$end),
                     strand=RRBS$strand, methylation = RRBS$methylation )

  }



  ####### only investigate bins that have RRBS cysosine methylation


  overlap.Granges.genomeVec = subsetByOverlaps( Granges.genomeVec, G.RRBS , minoverlap = 1L, ignore.strand = T  ) ; length(overlap.Granges.genomeVec)
  overlap.Granges.RRBS = subsetByOverlaps( G.RRBS ,  Granges.genomeVec, minoverlap = 1L, ignore.strand = T  ) ; length(overlap.Granges.RRBS)


  ###### Calculated each bin's average RRBS.CpG's methylation results as well as count how many 'hits' with RRBS.CpGs

  #### use IRanges objects for fast speed:

  I.overlap.genomeVec = IRanges(start=start(ranges(overlap.Granges.genomeVec)), end=end(ranges(overlap.Granges.genomeVec)))

  I.RRBS = IRanges(start=start(ranges(overlap.Granges.RRBS)), end=end(ranges(overlap.Granges.RRBS)))   # IRanges objects


  RRBS.methylation = overlap.Granges.RRBS$methylation

  overlap.bin.methyl = vector()
  num_bin_hit = vector()

  total.overlap.bin = length(overlap.Granges.genomeVec)

  cat( paste0("Calculating ", method, " correlation coefficient for chromosome ", CHR, ", bin size ", ws, " ... \n "  )      )

  for ( i in 1: total.overlap.bin  ) {

    hit = findOverlaps(  I.overlap.genomeVec[i,], I.RRBS, minoverlap = 1L   )

    num_bin_hit[i] = length(hit)  # number of RRBS.CpG hit by the bin

    overlap.bin.methyl[i] = mean( RRBS.methylation[ hit@to ] )
  }


  ### add to object overlap.Granges.genomeVec

  overlap.Granges.genomeVec$hits = num_bin_hit

  overlap.Granges.genomeVec$RRBS.methylation = overlap.bin.methyl



  #### Get adjusted methylation level from MeDIPS and compared genome correlation:

  ## Only use bins with at least 'minRRBS' number of RRBS CpGs for correlation calculation:


  hits.Granges.genomeVec = overlap.Granges.genomeVec[overlap.Granges.genomeVec$hits >= minRRBS, ]



  ####### result : correlation

  if ( method == "pearson" ) {

    cat( paste0( method, " correlation coefficient for chromosome ", CHR, ", using bin size ", ws, " is: \n", cor( hits.Granges.genomeVec$binMethyl , hits.Granges.genomeVec$RRBS.methylation, method = "pearson"),"\n" ) )

    correlation[as.character(CHR)] = cor( hits.Granges.genomeVec$binMethyl , hits.Granges.genomeVec$RRBS.methylation, method = "pearson")
  }

  if ( method == "spearman" ) {

    cat( paste0( method, " correlation coefficient for chromosome ", CHR, ", using bin size ", ws, " is: \n" ), cor( hits.Granges.genomeVec$binMethyl , hits.Granges.genomeVec$RRBS.methylation, method = "spearman"),"\n"  )

    correlation[as.character(CHR)] = cor( hits.Granges.genomeVec$binMethyl , hits.Granges.genomeVec$RRBS.methylation, method = "spearman")
  }

  if ( method == "kendall" ) {

    cat( paste0( method, " correlation coefficient for chromosome ", CHR, ", using bin size ", ws, " is: \n" ), cor( hits.Granges.genomeVec$binMethyl , hits.Granges.genomeVec$RRBS.methylation, method = "kendall"),"\n"  )

    correlation[as.character(CHR)] = cor( hits.Granges.genomeVec$binMethyl , hits.Granges.genomeVec$RRBS.methylation, method = "kendall")

  }

  }

  return(correlation)


}

