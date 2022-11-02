## UMI_barcodes
Full list of commonly used the UMI barcodes for cfMeDIP-seq
1) NNT_barcodes.txt           ## Barcodes for the pattern of NNT
2) UMI_barcodes_OICR.txt      ## Barcodes list applied by the OICR protocols

## Spike-in_genomes
1) SyntheticDNA_Arabidopsis_BACs.fa consists of Arabidopsis BAC (F19K16_F24B22) and [sythetic DNA sequences](https://github.com/hoffmangroup/2020spikein/tree/master/Preprocessing).

2) SyntheticDNA_Arabidopsis_BACs_seqNames.txt: sequences' name

### How to forge a BSgenome package for the spike-ins
1) Spike-in genome
  ```bash
  ## Get Fasta sequence and transfer to 2bit format with ucsctools
  $ faToTwoBit BCA_F19K16_F24B22.fa BCA_F19K16_F24B22.2bit
  ```
2) Forge BSgenome package
  ```{r }
  # prepare the seed file according to BSgenome instruction
  # eg: BSgenome.Athaliana.BAC.F19K16.F24B22-seed
  library(BSgenome)
  forgeBSgenomeDataPkg("path/to/seed/file")
  ```
3) Build package
  ```bash
  $ R CMD build /path/to/pkgdir
  ```  

## Filter_regions
Gaps and filters regions to be removed in fragment profile analysis!!  
1) gaps_filters_hg19.rdata and gaps_filters_hg38.rdata were produced by the script gaps_filters_hg.R, which includes telomeres, centromeres and ENOCDE blacklist regions!!
2) AB_hg19.rdata and ab_hg38.rdata are HiC_AB_Compartments downloaded & liftovered from [here](https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt)
