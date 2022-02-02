#!/bin/bash

## the script will build BWA index for combined human and spike-in genomes.
## "Usage: ./download_build_reference.sh [GENOME] [SPIKEIN_FA] [INDEX_PREFIX] [DEST_DIR]"
## "Example: ./download_build_reference.sh hg38 ./data/BAC_F19K16_F24B22.fa hg38_BAC_F19K16_F24B22 /your/genome/data/path/hg38"


#################
## initilizaiton
#################
source ~/miniconda3/etc/profile.d/conda.sh

# OR:
# source ~/conda3/etc/profile.d/conda.sh

TSV=${DEST_DIR}/${GENOME}.tsv

## human genome
hg_fa=` awk '{if($1=="ref_fa") print $2}' $TSV`


## combine genSPIKEIN_FAomes

cat ${hg_fa} ${SPIKEIN_FA} > ${DEST_DIR}/${hg38_BAC_F19K16_F24B22}.fa


#################################################
## build bwa index for combined reference genomes
#################################################
## combine genomes
cat $(basename ${REF_FA}) TAIR10_chr_all.fas.gz > ${GENOME}_tair10.fa.gz


## Might need to build index separatly due to the storage limitation
## in login homefolder (50G quota for H4H)

echo "=== Building bwa index for mereged genomes ..."
conda activate tcge-cfmedip-seq-pipeline

#bwa index -a bwtsw ${GENOME}_tair10.fa

#mkdir -p bwa_index_${GENOME}_tair10
#mv ${GENOME}_tair10* ./bwa_index_${GENOME}_tair10

## bwa index merged
#bwa_idx_merged=$(ls $PWD/bwa_index_${GENOME}_tair10/*fa)
#echo -e "bwa_idx_hg_tair\t${bwa_idx_merged}" >> ${TSV}

echo "=== Done! ==="
