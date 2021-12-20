##################################################
##### conda_env.yaml and conda_env_r.yaml  #######

## conda_env.ymal
conda install -c anaconda graphviz

conda env export -c bioconda  --from-history | grep -v "^prefix" >  conda_env.yaml


## conda_env_r.yaml

conda install -c bioconda bioconductor-medips

conda env export -c bioconda  --from-history | grep -v "^prefix" >  conda_env_r.yaml


## export environment


###############################
##### test run on H4H  #######

conda activate cfmedip-seq-pipeline

## for h4h
## need to excute test run in workdir to install conda_env_R.yaml in .snakmake
## on the build node with internet

## generate the DGA plots ####
## dot -Tsvg

################################################################################
## for paired-end inputs
snakemake --snakefile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/config/config_pe_template.yaml \
          -p --dag | dot -Tpdf > /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/figures/dag_pe.pdf

## for single-end inputs
snakemake --snakefile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/config/config_se_template.yaml \
          -p --dag | dot -Tpdf > /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/figures/dag_se.pdf


#### test run on H4H without sbatch
snakemake --snakefile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/config/config_pe_template.yaml \
          --use-conda -np

#### test run on H4H with sbatch
## using the --unlock flag to remove a wkdir lock
snakemake --snakefile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/config/config_pe_template.yaml \
          --cluster-config /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/config/cluster_std_err.json \
          --use-conda --core 2 --cluster "sbatch -p all -o {cluster.std} -e {cluster.err}" \
          --latency-wait 10 --jobs 8 -np





################################################################################
#### test run with cfmeDIP-seq data
snakemake --snakefile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/projects/tcge/cell_free_epigenomics/test_run/config_test.yaml \
          --use-conda -np


## run on cluster by sbatch  :: multiple jobs will submitted
snakemake --snakefile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/projects/tcge/cell_free_epigenomics/test_run/config_test.yaml \
          --use-conda --cores 2 --cluster "sbatch -p all" --jobs 8
