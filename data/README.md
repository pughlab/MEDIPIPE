## UMI barcode list
1) NNT_barcodes.txt           ## full list of barcodes for the pattern of NNT
2) UMI_barcodes_OICR.txt      ## Barcodes list applied by the OICR protocols

## forge the spike-in BSgenome packages
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

## Gaps and filters regions to be removed in fragment profile analysis  
1) gaps_filters_hg19.rdata and gaps_filters_hg38.rdata were produced by the script gaps_filters_hg.R
2) AB_hg19.rdata and ab_hg38.rdata were HiC_AB_Compartments downloaded & liftover from [here](https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt)
