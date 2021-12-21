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

## --jobs   ## number of samples
## --cores  ## jobs*2

snakemake --snakefile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/config/config_pe_template.yaml \
          --use-conda  --conda-prefix /cluster/home/yzeng/miniconda3/envs/cfmedip-seq-pipeline_R \
          --cluster-config /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/config/cluster_std_err.json \
          --cluster "sbatch -p all --mem=16G -o {cluster.std} -e {cluster.err}" \
          --latency-wait 60 --cores 8 --jobs 4 -p


## move all submission std and and err to logs
mv submit* ./logs

conda deactivate
