#!/bin/bash

## the script was modified based on :
## https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/scripts/download_genome_data.sh
## ENCODE DCC Version: v3 for >=ENCODE4

## "A TSV file [DEST_DIR]/[GENOME].tsv will be generated. Use it for pipeline."
## "Supported genomes: hg19 and hg38"; Arabidopsis TAIR10 genome will be downloaded,
##  as well as building bwa index for merged genomes.
## "Usage: ./assets/Reference/download_build_reference.sh [GENOME] [DEST_DIR]"
## "Example: ./assets/Reference/download_build_reference.sh hg38 /your/genome/data/path/hg38"


#################
## initilizaiton
#################
source ~/miniconda3/etc/profile.d/conda.sh

# OR:
# source ~/conda3/etc/profile.d/conda.sh

GENOME=$1
DEST_DIR=$(cd $(dirname $2) && pwd -P)/$(basename $2)
TSV=${DEST_DIR}/${GENOME}.tsv

echo "=== Creating destination directory and TSV file..."
mkdir -p ${DEST_DIR}
cd ${DEST_DIR}

############################
## URL for reference genomes
############################

if [[ "${GENOME}" == "hg38" ]]; then
  REGEX_BFILT_PEAK_CHR_NAME="chr[\dXY]+"
  MITO_CHR_NAME="chrM"
  REF_FA="https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"
  REF_MITO_FA="https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15_mito_only/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15_mito_only.fasta.gz"
  BWT2_IDX="https://www.encodeproject.org/files/ENCFF110MCL/@@download/ENCFF110MCL.tar.gz"
  BWT2_MITO_IDX="https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15_mito_only_bowtie2_index/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15_mito_only_bowtie2_index.tar.gz"
  BWA_IDX="https://www.encodeproject.org/files/ENCFF643CGH/@@download/ENCFF643CGH.tar.gz"
  BWA_MITO_IDX="https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15_mito_only_bwa_index/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15_mito_only_bwa_index.tar.gz"
  CHRSZ="https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv"
  BLACKLIST="https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
  TSS="https://www.encodeproject.org/files/ENCFF493CCB/@@download/ENCFF493CCB.bed.gz"
  DNASE="https://www.encodeproject.org/files/ENCFF304XEX/@@download/ENCFF304XEX.bed.gz"
  PROM="https://www.encodeproject.org/files/ENCFF140XLU/@@download/ENCFF140XLU.bed.gz"
  ENH="https://www.encodeproject.org/files/ENCFF212UAV/@@download/ENCFF212UAV.bed.gz"

  #REF_FA_TAIR10="https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas"

fi

if [[ "${GENOME}" == "hg19" ]]; then
  REGEX_BFILT_PEAK_CHR_NAME="chr[\dXY]+"
  MITO_CHR_NAME="chrM"
  REF_FA="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz"
  REF_MITO_FA="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/male.hg19.chrM.fa.gz"
  CHRSZ="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/hg19.chrom.sizes"
  BWT2_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/bowtie2_index/male.hg19.fa.tar"
  BWT2_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/bowtie2_index/male.hg19.chrM.fa.tar"
  BWA_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/bwa_index/male.hg19.fa.tar"
  BWA_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/bwa_index/male.hg19.chrM.fa.tar"
  BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
  TSS="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/hg19_gencode_tss_unique.bed.gz"
  DNASE="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
  PROM="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
  ENH="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/reg2map_honeybadger2_dnase_enh_p2.bed.gz"

  ## Arabidopsis
  # REF_FA_TAIR10="https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas"

fi



####################
## Downloading files
####################

echo "=== Downloading files..."
wget -c -O $(basename ${REF_FA}) ${REF_FA}
wget -c -O $(basename ${REF_MITO_FA}) ${REF_MITO_FA}
wget -c -O $(basename ${CHRSZ}) ${CHRSZ}

## TAIR10
#wget -c -O $(basename ${REF_FA_TAIR10}) ${REF_FA_TAIR10}
#sed -i -e 's/^>/>tair10_chr/' TAIR10_chr_all.fas
#gzip  TAIR10_chr_all.fas

## combine genomes
# cat $(basename ${REF_FA}) TAIR10_chr_all.fas.gz > ${GENOME}_tair10.fa.gz

## annotated regions
wget -N -c ${BLACKLIST}
wget -N -c ${TSS}
wget -N -c ${DNASE}
wget -N -c ${PROM}
wget -N -c ${ENH}

###  to add CpG islands, sea, shore ..

gunzip *.gz

echo "=== Downloading bwa index..."

if [[ "${GENOME}" == "hg38" ]]; then
mkdir -p ${DEST_DIR}/bwa_index_${GENOME}
cd ${DEST_DIR}/bwa_index_${GENOME}
wget -c ${BWA_IDX}
tar xvzf *.tar.gz
rm *.tar.gz
fi

if [[ "${GENOME}" == "hg19" ]]; then
mkdir -p ${DEST_DIR}/bwa_index_${GENOME}
cd ${DEST_DIR}/bwa_index_${GENOME}
wget -c ${BWA_IDX}
tar xf *.tar
#rm *.tar
fi

#################################################
## build bwa index for combined reference genomes
#################################################

## Might need to build index separatly due to the storage limitation
## in login homefolder (50G quota for H4H)
## echo "=== Building bwa index for mereged genomes ..."

#conda activate tcge-cfmedip-seq-pipeline
#bwa index -a bwtsw ${GENOME}_tair10.fa

#mkdir -p bwa_index_${GENOME}_tair10
#mv ${GENOME}_tair10* ./bwa_index_${GENOME}_tair10


####################
## Creating TSV file
####################
cd ${DEST_DIR}
rm -f ${TSV}
touch ${TSV}

echo -e "genome_name\t${GENOME}" >> ${TSV}
echo -e "ref_fa\t${DEST_DIR}/$(basename ${REF_FA})" >> ${TSV}
echo -e "ref_mito_fa\t${DEST_DIR}/$(basename ${REF_MITO_FA})" >> ${TSV}
echo -e "chrsz\t${DEST_DIR}/$(basename ${CHRSZ})" >> ${TSV}
echo -e "blacklist\t${DEST_DIR}/$(basename ${BLACKLIST})" >> ${TSV}
echo -e "tss\t${DEST_DIR}/$(basename ${TSS})" >> ${TSV}
echo -e "dnase\t${DEST_DIR}/$(basename ${DNASE})" >> ${TSV}
echo -e "prom\t${DEST_DIR}/$(basename ${PROM})" >> ${TSV}
echo -e "enh\t${DEST_DIR}/$(basename ${ENH})" >> ${TSV}
#echo -e "ref_fa_tair10\t${DEST_DIR}/$(basename ${REF_FA_TAIR10})" >> ${TSV}

## bwa index
if [[ "${GENOME}" == "hg38" ]]; then
bwa_idx=$(ls $PWD/bwa_index_${GENOME}/*fna)
echo -e "bwa_index\t${bwa_idx}" >> ${TSV}
fi

if [[ "${GENOME}" == "hg19" ]]; then
mv male.hg19.fa ./bwa_index_${GENOME}
bwa_idx=$(ls $PWD/bwa_index_${GENOME}/*fa)
echo -e "bwa_index\t${bwa_idx}" >> ${TSV}
fi

## bwa index merged
#bwa_idx_merged=$(ls $PWD/bwa_index_${GENOME}_tair10/*fa)
#echo -e "bwa_idx_spikein\t${bwa_idx_merged}" >> ${TSV}

## remove gz suffix
sed -i 's/.gz//g' ${TSV}

echo "=== Done! ==="
