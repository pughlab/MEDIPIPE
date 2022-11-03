#!/bin/bash
#SBATCH -p all              ## Specify SLURM partition for job submssion
#SBATCH -t 5-00:00:00
#SBATCH --mem=16G
#SBATCH -J submit_snakemake_%j
#SBATCH -o submit_snakemake_%j.out
#SBATCH -e submit_snakemake_%j.err

## submit the snakemake for sample list

## configure shell: full path to conda.sh
source ~/miniconda3/etc/profile.d/conda.sh

conda activate MEDIPIPE

## mkdir for cluster submission logs
## defined in .workflow/config/cluster_std_err.json
cd  /path/to/test/Res
mkdir -p logs_cluster

## unlock workdir just in case the folder locked accidently before
snakemake --snakefile /path/to/workflow/Snakefile \
          --configfile /path/to/test/config_template.yaml \
          --unlock

## --jobs   ## number of samples
## sbatch -c :  maximal 8 threads per multithreading job by default, less -c INT  will be scaled down to INT

snakemake --snakefile /path/to/workflow/Snakefile \
          --configfile /path/to/test/config_template.yaml \
          --use-conda  --conda-prefix ${CONDA_PREFIX}_extra_env \
          --cluster-config /path/to/workflow/config/cluster_std_err.json \
          --cluster "sbatch -p long -c 8 --mem=16G -o {cluster.std} -e {cluster.err}" \
          --latency-wait 60 --jobs 4 -p

conda deactivate
