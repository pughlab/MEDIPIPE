## UMI barcode list


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
