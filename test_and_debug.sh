##################################################
##### conda_env.yaml and conda_env_R.yaml  #######

## conda_env.ymal
conda install -c anaconda graphviz
...
conda env export -c bioconda  --from-history | grep -v "^prefix" >  conda_env.yaml


###############################
##### test run on H4H  #######
###############################
conda activate tcge-cfmedip-seq-pipeline


########################################
## generate the DGA plots with dot -Tsvg
## for paired-end inputs
snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/config_template.yaml \
          -p --dag | dot -Tpdf > /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/figures/dag.pdf     ## dag_pe | dag_se



############################
#### test run with test data
############################

##  stand-alone dry-run
snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/config_template.yaml \
          --unlock

snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/config_template.yaml \
          --use-conda  --conda-prefix ~/miniconda3/envs/tcge-cfmedip-seq-pipeline-sub \
          --cores 4 -pn

## sbatch /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/sbatch_snakemake_template.sh
## using the --unlock flag to unlock workdir
snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/config_template.yaml \
          --use-conda  --conda-prefix /cluster/home/yzeng/miniconda3/envs/tcge-cfmedip-seq-pipeline-sub \
          --cluster-config /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/cluster_std_err.json \
          --cluster "sbatch -p all --mem=16G -o {cluster.std} -e {cluster.err}" \
          --latency-wait 60 --cores 12 --jobs 4 -p

#############################################
#### test run with real tcge-cfmedip-seq data
#############################################
snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/config_real_run.yaml \
          --unlock

##  stand-alone dry-run
snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/config_real_run.yaml \
          --use-conda  --conda-prefix /cluster/home/yzeng/miniconda3/envs/tcge-cfmedip-seq-pipeline-sub \
          --cores 4 -pn

snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/config_real_run.yaml \
          -p --dag | dot -Tpdf > /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/figures/real_run.pdf

## sbatch /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/sbatch_snakemake_real_run.sh
snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/config_real_run.yaml \
          --use-conda  --conda-prefix /cluster/home/yzeng/miniconda3/envs/tcge-cfmedip-seq-pipeline-sub \
          --cluster-config /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/config/cluster_std_err.json \
          --cluster "sbatch -p himem -c 12 --mem=60G -J {cluster.jid} -o {cluster.std} -e {cluster.err}" \
          --latency-wait 60 --jobs 8 -p


###############################
##### testing single end ######
###############################
snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/projects/tcge/cell_free_epigenomics/processed_data/config_MCA.yaml \
          --use-conda  --conda-prefix /cluster/home/yzeng/miniconda3/envs/tcge-cfmedip-seq-pipeline-sub \
          --cores 4 -pn

#######################################
## testing sample with  UMIs
snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/projects/tcge/cell_free_epigenomics/processed_data/config_HNSC.yaml \
          --unlock

snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/projects/tcge/cell_free_epigenomics/processed_data/config_HNSC.yaml \
          --use-conda  --conda-prefix /cluster/home/yzeng/miniconda3/envs/tcge-cfmedip-seq-pipeline-sub \
          --cores 4 -pn


###########################################
##### testing samples with spike ins + UMI  ######
#####  and UMIs  data
###########################################
snakemake --snakefile /cluster/home/yzeng/snakemake/tcge-cfmedip-seq-pipeline/workflow/Snakefile \
          --configfile /cluster/projects/tcge/cell_free_epigenomics/processed_data/config_PDAC.yaml \
          --use-conda  --conda-prefix /cluster/home/yzeng/miniconda3/envs/tcge-cfmedip-seq-pipeline-sub \
          --cores 4 -pn
