#!/bin/bash

## the script will build BWA index for combined human and spike-in genomes.
## "Usage: ./build_reference_index.sh [GENOME] [SPIKEIN_FA] [INDEX_PREFIX] [DEST_DIR]"
## "Example: ./build_reference_index.sh hg38 ./data/BAC_F19K16_F24B22.fa hg38_BAC_F19K16_F24B22 /your/genome/data/path/hg38"

## "Example: ./build_reference_index.sh hg38 ./data/BAC_F19K16_F24B22.fa hg38_BAC_F19K16_F24B22 /cluster/projects/tcge/DB/cfmedip-seq-pepeline/hg38"

#################
## initilizaiton
#################
source ~/miniconda3/etc/profile.d/conda.sh

# OR:
# source ~/conda3/etc/profile.d/conda.sh

GENOME=$1
SPIKEIN_FA=$2
INDEX_PREFIX=$3
DEST_DIR=$4

## human genome
TSV=${DEST_DIR}/${GENOME}.tsv
hg_fa=` awk '{if($1=="ref_fa") print $2}' $TSV`

## combine human genome with sipke-in sequence
cat ${hg_fa} ${SPIKEIN_FA} > ${DEST_DIR}/${INDEX_PREFIX}.fa


#################################################
## build bwa index for combined reference genomes
#################################################


## Might need to build index separatly due to the storage limitation
cd ${DEST_DIR}

echo "=== Building bwa index for mereged genomes ..."
conda activate tcge-cfmedip-seq-pipeline

bwa index -a bwtsw ${INDEX_PREFIX}.fa

mkdir -p bwa_index_${INDEX_PREFIX}
mv ${INDEX_PREFIX}.* ./bwa_index_${INDEX_PREFIX}

## adding new bwa index into manifest file
bwa_idx_merged=$(ls $PWD/bwa_index_${INDEX_PREFIX}/*fa)
echo -e "bwa_idx_hg_spikein\t${bwa_idx_merged}" >> ${TSV}

echo "=== Done! ==="
