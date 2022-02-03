
# Dependencies

The R/Pythone packages are unavailable through conda/bioconda/r were cloned.

1. [ConsensusCruncher](https://github.com/pughlab/ConsensusCruncher) was cloned here to deal with the sequencing libraries come with custom barcodes, for instance, the Unique Molecular Identifier (UMI) enables accurate bioinformatic identification of PCR duplicates. A text file with these barcode sequences, one per line, is required for ConsensusCruncher's `--blist` parameter.  

2. [MeDEStrand] (https://github.com/eyzhao/MeDEStrand) was built on top of [MEDIPS](https://bioconductor.riken.jp/packages/3.8/bioc/html/MEDIPS.html), enabling estimation of absolute DNA methylation levels after correcting the CpG bias from the enrichment results.
