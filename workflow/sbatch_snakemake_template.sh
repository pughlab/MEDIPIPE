#!/bin/bash
#SBATCH -p all              ## for compile jobs and submit them
#SBATCH -t 5-00:00:00
#SBATCH --mem=16G
#SBATCH -J submit_snakemake_%j
#SBATCH -o submit_snakemake_%j.out
#SBATCH -e submit_snakemake_%j.err

## submit the snakemake for sample list
## sbatch /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/sbatch_snakemake_template.sh

# srun snakemake
conda activate cfmedip-seq-pipeline

snakemake --snakefile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/config/config_pe_template.yaml \
          --cluster-config /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/config/cluster_std_err.json \
          --cluster "sbatch -p all --mem=16G -o {cluster.std} -e {cluster.err}" \
          --latency-wait 10 --use-conda --core 2 --jobs 8 -p

## move all submission std and and err to logs
mv submit* ./logs

conda deactivate
