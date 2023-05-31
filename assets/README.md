## Filter_regions
Gaps and filters regions to be removed in fragment profile analysis!!  
1) gaps_filters_hg19.rdata and gaps_filters_hg38.rdata were produced by the script gaps_filters_hg.R, which includes telomeres, centromeres and ENOCDE blacklist regions!!
2) AB_hg19.rdata and ab_hg38.rdata are HiC_AB_Compartments downloaded & liftovered from [here](https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt)

## Reference
* Download reference genome data
You can download reference genome, pre-build BWA index and annotated regions (e.g., blacklist) from ENCODE for hg38 and hg19 on the command line. The manifest file hg38/hg19.tsv will be generated accordingly. Currently, the ENCODE black list and bwa index are mandatory for the manifest file, which you can also create it by yourself based on `.Reference/hg38_template.tsv` with existing data.


```bash
## eg: ./download_build_reference.sh hg38 /your/genome/data/path/hg38
$ ./assets/Reference/download_reference.sh [GENOME] [DEST_DIR]
```

* Build reference genomes index
If your sequencing libraries come with spike-ins, you can build new aligner index after combining spike-in genome with human genome. The new index information will be appended to corresponding manifest file.

```bash
## eg: ./assets/Reference/build_reference_index.sh hg38 ./data/BAC_F19K16_F24B22.fa hg38_BAC_F19K16_F24B22 /your/genome/data/path/hg38
$ ./assets/Reference/build_reference_index.sh [GENOME] [SPIKEIN_FA] [INDEX_PREFIX] [DEST_DIR]
```


## Spike-in_genomes
Spike-in FASTA sequences for two BACs: [F19K16](https://www.arabidopsis.org/servlets/TairObject?type=assembly_unit&id=362) from Arabidopsis Chr1 and [F24B22](https://www.arabidopsis.org/servlet/TairObject?type=AssemblyUnit&name=F24B22) from Arabidopsis Chr3, and [sytheticDNAs](https://github.com/hoffmangroup/2020spikein) were enclosed.

* SyntheticDNA_Arabidopsis_BACs.fa consists of Arabidopsis BAC (F19K16_F24B22) and [sythetic DNA sequences](https://github.com/hoffmangroup/2020spikein/tree/master/Preprocessing).

* SyntheticDNA_Arabidopsis_BACs_seqNames.txt: sequences' name

* How to forge a BSgenome package for the spike-ins
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

## UMI_barcodes
Full list of commonly used the UMI barcodes for cfMeDIP-seq
1) NNT_barcodes.txt           ## Barcodes for the pattern of NNT
2) UMI_barcodes_OICR.txt      ## Barcodes list applied by the OICR protocols  
