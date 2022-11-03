## Input configfile specification
> **IMPORTANT**: READ THROUGH THE GUIDE INFORMATION IN THE TEMPLATE TO MAKE A CORRECT CONFIG FILE. ESPECIALLY FOR SAMPLES WITH SPIKE-IN CONTROL AND/OR UMI BARCODES.

A config YAML file specifies all PATHs of input files and parameters that are necessary for successfully running this pipeline. This includes a specification of the path to the genome reference files. Please make sure to specify absolute paths rather than relative paths in your config files. More detail can be found at [here](./test/config_template.yaml)

1) The samples's sequence data table template. Note:Prepare the table for single-end and paired-end samples separately and use exactly same table `header`, if there are multiple lanes, use comma to separate the list.

|	sample_id   |     R1	     |  R2(p.r.n.)|
|-------------|--------------|------------|
|  A	|  full/path/to/A_L001_R1.fq.gz |                              |
|  B	|  full/path/to/B_L001_R1.fq.gz | full/path/to/B_L001_R2.fq.gz |
|  C  |  path/C_L001_R1.fq.gz,path/C_L002_R1.fq.gz | path/C_L001_R2.fq.gz,path/C_L002_R2.fq.gz  |






## test-datasets: `atacseq`

This test data (ATAC-seq) was download from [nf-core/atacseq](https://github.com/nf-core/test-datasets/tree/0c58a9f36205cc5f8c6bbb2ca03c401c61cb849d).

`Reference/`: Genome reference and BWA index (iGenomes sacCer3 release)   
`Fastq/` : FastQ files sub-sampled to 100,000 paired-end reads   

## test dataset origin:

*S. cerevisiae* paired-end ATAC-seq dataset was obtained from:

Schep AN, Buenrostro JD, Denny SK, Schwartz K, Sherlock G, Greenleaf WJ. Structured nucleosome fingerprints enable high-resolution mapping of chromatin architecture within regulatory regions. Genome Res 2015 Nov;25(11):1757-70. [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/26314830) [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66386)

### Sample information

| GEO_ID	    |   LAYOUT	  | SAMPLE_NAME|
|-------------|-------------|------------|
| SRR1822153	| paired-end	|     A	     |
| SRR1822154	| paired-end	|     B	     |
| SRR1822157	| single-end	|     C      |
| SRR1822158	| single-end	|     D	     |
