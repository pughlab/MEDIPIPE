#!/bin/bash
#SBATCH -p all              ## for compile jobs and submit them
#SBATCH -t 5-00:00:00
#SBATCH --mem=16G
#SBATCH -J submit_snakemake_%j
#SBATCH -o submit_snakemake_%j.out
#SBATCH -e submit_snakemake_%j.err

## submit the snakemake for sample list
## sbatch /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/sbatch_snakemake_template.sh

## configure shell: full path to conda.sh
source ~/miniconda3/etc/profile.d/conda.sh

conda activate tcge-cfmedip-seq-pipeline

## mkdir for cluster submission logs
## defined in .workflow/config/cluster_std_err.json
cd  /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline-test-run
mkdir -p logs_cluster

## unlock workdir just in case the folder locked accidently before
snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/config_template.yaml \
          --unlock

## --jobs   ## number of samples
## sbatch -c :  maximal 8 threads per multithreading job by default, less -c INT  will be scaled down to INT

snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/config_template.yaml \
          --use-conda  --conda-prefix /cluster/home/yzeng/miniconda3/envs/tcge-cfmedip-seq-pipeline-sub \
          --cluster-config /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/cluster_std_err.json \
          --cluster "sbatch -p long -c 8 --mem=16G -o {cluster.std} -e {cluster.err}" \
          --latency-wait 60 --jobs 4 -p

conda deactivate
