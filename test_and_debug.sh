##################################################
##### conda_env.yaml and conda_env_r.yaml  #######

## conda_env.ymal


conda install -c anaconda graphviz

conda env export -c bioconda  --from-history | grep -v "^prefix" >  conda_env.yaml


## conda_env_r.yaml

conda install -c bioconda bioconductor-medips

conda env export -c bioconda  --from-history | grep -v "^prefix" >  conda_env_r.yaml


## export environment


conda activate cfmedip-seq-pipeline



###############################
##### test run on H4H  #######

conda activate cfmedip-seq-pipeline
snakemake --snakefile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile --create-envs-only --use-conda

# run test run to install R env as well tesing working

#testing using conda-prefix
snakemake --snakefile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/config/config_pe_template.yaml \
          --use-conda  --conda-prefix /cluster/home/yzeng/miniconda3/envs/cfmedip-seq-pipeline_R -p


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
          --use-conda  --conda-prefix /cluster/home/yzeng/miniconda3/envs/cfmedip-seq-pipeline_R \
          --cores 4 -p


#### test run on H4H with sbatch
## using the --unlock flag to remove a wkdir lock
snakemake --snakefile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/config/config_pe_template.yaml \
          --use-conda  --conda-prefix /cluster/home/yzeng/miniconda3/envs/cfmedip-seq-pipeline_R \
          --cluster-config /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/config/cluster_std_err.json \
          --cluster "sbatch -p all --mem=16G -o {cluster.std} -e {cluster.err}" \
          --latency-wait 60 --cores 8 --jobs 4 -np







################################################################################
#### test run with cfmeDIP-seq data
snakemake --snakefile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/projects/tcge/cell_free_epigenomics/test_run/config_test.yaml \
          --use-conda -np


## run on cluster by sbatch  :: multiple jobs will submitted
snakemake --snakefile /cluster/home/yzeng/snakemake/cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/projects/tcge/cell_free_epigenomics/test_run/config_test.yaml \
          --use-conda --cores 2 --cluster "sbatch -p all" --jobs 8
