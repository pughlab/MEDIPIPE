# MEDIPIPE: (cf)MeDIP-seq Data QC and Analysis Pipeline

## Intoduction(v1.0.0)

The MEDIPIPE is designed for automated end-to-end quality control and processing of cfMeDIP-seq and MeDIP-seq data. The pipeline was developed with [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html), which will automatically deploy the execution environments. The pipeline can be run on compute clusters with job submission engines as well as on stand along machines. The pipeline starts from the raw FASTQ files all the way to quality metrics and methylation quantification/estimation matrixes. The pipeline supports both signle-end and paired-end sequencing data with or without spike-in/UMI sequences. The pipeline can be applied to individual sample, as well as to aggregate multiple samples' to create matrixes of QC metrics; fragment profiles and methylation quantifications.

The pipeline was developed by [Yong Zeng](mailto:yzeng@uhnresearch.ca) based on some prior work of Wenbin Ye, [Eric Zhao](https://github.com/pughlab/cfMeDIP-seq-analysis-pipeline). The pipeline is available on [COBE](https://www.pmcobe.ca/pipeline/624d0bb4002b11003426f7d8) as well.


### Features

- **Portability**: The pipeline run can be performed across different cluster engines such as SLURM (tested). For other platforms, please refer to [here](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).
- **Supported genomes**: We provide genome database, which includes aligner indices and black list, downloader for human hg38, hg19 and along with Arabidopsis thaliana genome TAIR10. In addition, spike-in fasta sequences for two BACs: [F19K16](https://www.arabidopsis.org/servlets/TairObject?type=assembly_unit&id=362) from Arabidopsis Chr1 and [F24B22](https://www.arabidopsis.org/servlet/TairObject?type=AssemblyUnit&name=F24B22) from Arabidopsis Chr3, and [sytheticDNAs](https://github.com/hoffmangroup/2020spikein) were enclosed in `data/SyntheticDNA_Arabidopsis_BACs.fa`. You can also build genome database from FASTA for your custom genomes.


### How it works
This schematic diagram shows you how pipeline works:

<img src="figures/MEDIPIPE_Flowchart.png" alt="Schematic_diagram" style="width:100.0%" />

## Installation

1) Make sure that you have a Conda-based Python3 distribution(e.g.,the [Miniconda](https://docs.conda.io/en/latest/miniconda.html)). The installation of [Mamba](https://github.com/mamba-org/mamba) is also recommended since it is better at handling environments and installations:

	```bash
	$ conda install -n base -c conda-forge mamba
	```

2) Git clone this pipeline.
	```bash
	$ cd
	$ git clone https://github.com/yzeng-lol/MEDIPIPE
	```

3) Install pipeline\'s core enviroment
	```bash
	$ cd MEDIPIPE
	$ conda activate base
	$ mamba env create --file conda_env.yaml
	```

4) Test run
	> **IMPORTANT**: EXTRA ENVIRONMENTS WILL BE INSTALLED, MAKE SURE YOU STILL HAVE INTERNET ACCESS.
	* **Step 1:** Prepare reference, samples FASTQ and aggregation files according to templates in ./test.
	* **Step 2:** Specify input configuration file by following the instructions in ./test/config_template.yaml.
	* **NOTE:** For testing run, you can simply run the SED command below to specify files in Step1,2. The outputs can be found in ./test/Res.

	```bash
    $ sed -i 's,/path/to,'"$PWD"',g' ./test/\*template.\*
	$ conda activate MEDIPIPE
	$ snakemake --snakefile ./workflow/Snakefile \
	            --configfile ./test/config_template.yaml \
		    --conda-prefix ${CONDA_PREFIX}_extra_env \
	            --use-conda --cores 4 -pn
	```

5) Run on HPCs

	You can submit this pipeline on clusters after editing ./workflow/sbatch_snakemake_template.sh according different resource management system. The file must have the appropriate path names and parameters given the input data. The example here is for submitting a job with SLURM, however, this template could be modified according to different resource management systems. More details about cluster configuration can be found at [here](https://snakemake.readthedocs.io/en/stable/executing/cluster.html). For example:

	```bash
	## template is based on SLURM
	$ vim ./workflow/sbatch_snakemake_template.sh
	$ sbatch ./workflow/sbatch_snakemake_template.sh
	```

## Input configfile specification
> **IMPORTANT**: READ THROUGH THE GUIDE INFORMATION IN THE TEMPLATE TO MAKE A CORRECT CONFIG FILE. ESPECIALLY FOR SAMPLES WITH SPIKE-IN CONTROL AND/OR UMI BARCODES.

A config YAML file specifies all PATHs of input files and parameters that are necessary for successfully running this pipeline. This includes a specification of the path to the genome reference files. Please make sure to specify absolute paths rather than relative paths in your config files. More detail can be found at [here](./test/config_template.yaml)

## Assets
There are serveral several scripts are enclosed in the [Assets](./Assets/config_template.yaml), allowing you to download/build reference index and manifest table, to fogre BSgeome package for spike-in controls
