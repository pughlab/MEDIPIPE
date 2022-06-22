#######################################
## forge the spike-in BSgenome packages

## UCSCtools 
$ faToTwoBit BCA_F19K16_F24B22.fa BCA_F19K16_F24B22.2bit

# prepare the seed file accordingly
# BSgenome.Athaliana.BAC.F19K16.F24B22-seed

## in R 
> library(BSgenome)
> forgeBSgenomeDataPkg("path/to/seed/file")

## build package
R CMD build /path/to/pkgdir 
##########################################
