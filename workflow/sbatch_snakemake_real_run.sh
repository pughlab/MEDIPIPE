#!/bin/bash
#SBATCH -p all              ## for compile jobs and submit them
#SBATCH -t 5-00:00:00
#SBATCH --mem=16G
#SBATCH -J submit_snakemake_%j
#SBATCH -o submit_snakemake_%j.out
#SBATCH -e submit_snakemake_%j.err

## submit the snakemake for sample list
## sbatch /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/sbatch_snakemake_real_run.sh

# run snakemake
#conda init bash
conda activate tcge-cfmedip-seq-pipeline


## unlock workdir just in case the folder locked accidently before
snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/projects/tcge/cell_free_epigenomics/test_run/config_test.yaml \
          --unlock


## --jobs   ## number of samples
## --cores  ## jobs*2

snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/projects/tcge/cell_free_epigenomics/test_run/config_test.yaml \
          --use-conda  --conda-prefix /cluster/home/yzeng/miniconda3/envs/tcge-cfmedip-seq-pipeline-sub \
          --cluster-config /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/cluster_std_err.json \
          --cluster "sbatch -p all --mem=16G -o {cluster.std} -e {cluster.err}" \
          --latency-wait 500 --cores 8 --jobs 4 -p

## move all submission std and and err to logs
mv submit* ./logs

conda deactivate
