# TCGE cfMeDIP-seq pipeline

## Intoduction

The Cancer Genetics and Epigenetics (TCGE) program cfMeDIP-seq pipeline is designed for automated end-to-end quality control and processing of cfMeDIP-seq and MeDIP-seq data. The pipeline was developed with [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html), which will automatically deploy the execution environment. The pipeline can be run on compute clusters with job submission engines as well as on stand along machines. The pipeline starts from the raw FASTQ files all the way to QC matrix and signal track generation. The pipeline supports both signle-end and paired-end sequencing data with or without spike-in/UMI sequences. The outputs produced by the pipeline include 1) formatted QC reports, 2) enrichment signal tracks, 3)...

The pipeline was developed by Yong Zeng based on some prior works of Wenbin Ye, Eric Zhao, ..

### Features

- **Portability**: The pipeline run can be performed across different cluster engines such as SLURM,...
- **Supported genomes**: We provide genome database, which includes aligner indices and black list, downloader for human hg38, hg19 and along with Arabidopsis thaliana genome TAIR10. In addition, fasta sequence for two BACs: [F19K16](https://www.arabidopsis.org/servlets/TairObject?type=assembly_unit&id=362) Arabidopsis Chr1 and [F24B22](https://www.arabidopsis.org/servlet/TairObject?type=AssemblyUnit&name=F24B22) from Arabidopsis Chr3 was enclosed in `data/BAC_F19K16_F24B22.fa`. You can also build genome database from FASTA for your custom genomes.


### How it works
This schematic diagram shows you how pipeline works:

<img src="figures/cfmedip-seq-pipeline.png" alt="Schematic_diagram" style="width:100.0%" />

## Installation

1) Make sure that you have a Conda-based Python3 distribution. The recommended choice is [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge). In case you don't use Mambafore, you can always install [Mamba](https://github.com/mamba-org/mamba) into any other Conda-based Python distribution with:

	```bash
	$ conda install -n base -c conda-forge mamba
	```

2) Git clone this pipeline.
	```bash
	$ cd
	## token request for private repository
	$ git clone https://github.com/yzeng-lol/tcge-cfmedip-seq-pipeline
	```

3) Install pipeline\'s core enviroments
	```bash
	$ cd tcge-cfmedip-seq-pipeline
	$ conda activate base
	$ mamba env create --file conda_env.yaml
	```

4) Test run
	> **IMPORTANT**: sub envs will also be created during the test run.

	```bash
	$ conda activate tcge-cfmedip-seq-pipeline

	## then edit ./workflow/config/config_pe_template.yaml accordingly
	$ mkdir ../tcge-cfmedip-seq-pipeline-test-run
	$ vim ./workflow/config/config_test_run.yaml

	## run with the internet connection as well
	$ snakemake --snakefile ./workflow/Snakefile \
	            --configfile ./workflow/config/config_template.yaml \
		    --conda-prefix /path/to/conda/envs/tcge-cfmedip-seq-pipeline-sub \
	            --use-conda --cores 4 -p
	```


## Input files specification

### Download reference genome data
You can download reference genome, pre-build BWA index and annotated regions from ENCODE for hg38 and hg19 via following command line. The manifest file hg38/hg19.tsv will be generated accordingly.

```bash
## eg: ./download_build_reference.sh hg38 /your/genome/data/path/hg38
$ ./download_reference.sh [GENOME] [DEST_DIR]
```

### Build reference genomes index
If your sequencing libraries come with spike-ins, you can build new aligner index after combining spike-in genome with human genome. The new index information will be appended to corresponding manifest file.

```bash
## eg: ./build_reference_index.sh hg38 ./data/BAC_F19K16_F24B22.fa hg38_BAC_F19K16_F24B22 /your/genome/data/path/hg38
$ ./build_reference_index.sh [GENOME] [SPIKEIN_FA] [INDEX_PREFIX] [DEST_DIR]
```

### config file for input samples


## Run on HPCs

You can submit this pipeline on clusters after editing ./workflow/sbatch_snakemake_template.sh according different resource management system. More details about cluster configuration can be found at [here](https://snakemake.readthedocs.io/en/stable/executing/cluster.html). For example:

```bash
$ vim ./workflow/sbatch_snakemake_template.sh
$ sbatch ./workflow/sbatch_snakemake_template.sh
```
