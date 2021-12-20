#!/bin/bash
#SBATCH -p all              ## for compile jobs and submit them
#SBATCH -t 5-00:00:00
#SBATCH --mem=16G
#SBATCH -J sbatch_snakemake_%j
#SBATCH -o sbatch_snakemake_%j.out
#SBATCH -e sbatch_snakemake_%j.err


# srun snakemake
# --use-conda \
conda activate cfmedip-seq-pipeline

snakemake  -j 10 --use-conda \
  --snakefile  /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile \
  --configfile /cluster/projects/tcge/cell_free_epigenomics/test_run/config_test.yaml \
  --cluster-config /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/config/cluster_std_err.json \
  --cluster "sbatch -p all -o {cluster.std} -e {cluster.err}" -np

conda deactivate
