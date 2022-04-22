###############################
## setting up working directory
workdir: config['workdir']

## read in sample list
import pandas as pd

## read in sample and corresponding fq files talbe
SAMPLES = (
    pd.read_csv(config["samples"], sep="\t")
    .set_index("sample_id", drop=False)
    .sort_index()
)

## read in refrence files' info
REF = pd.read_csv(config["ref_files"], sep="\t", header = None, index_col = 0)
blacklist = REF.loc["blacklist"][1]   ## ENCODE blacklist


## paths for pipeline and/or reference data
wd = config["workdir"]
pipe_dir = config["pipeline_dir"]
umi_list = config["umi_list"]


#############################################
## get taget outputs based on the config file
#############################################

def get_rule_all_input():
    ## outputs for each individual sample
    #map_out = expand("sorted_reads/{sample}_sorted.bam.bai", sample = SAMPLES["sample"]),
    #dedup_out = expand("dedup_bam/{sample}_dedup.bam.bai", sample = SAMPLES["sample_id"]),
    #meth_qc = expand("meth_qc_quant/{sample}_meth_qc.txt", sample = SAMPLES["sample_id"]),
    #meth_quant = expand("meth_qc_quant/{sample}_meth_quant.RData", sample = SAMPLES["sample_id"]),

    ## aggregated outputs for SAMPLES
    meth_qc = "aggregated/meth_qc.txt",
    meta_quant = "aggregated/meth_count.txt.gz",
    meth_filt = "autos_bfilt/meth_count_autos_bfilt.txt.gz",

    ## single-end or paired-end
    if config["paired-end"]:
        ## FASQC out for raw and trimmed paired-end fqs (individual sample)
        #fastqc_raw_pe = expand("fastqc/pe/{sample}_{mate}_fastqc.zip",
        #                        sample = SAMPLES["sample_id"],
        #                        mate = ["R1", "R2"]),
        mult_qc = "aggregated/QC_pe/multiqc_report.html",
        return  mult_qc + meth_qc + meta_quant + meth_filt

    else:
        ## FASQC out for raw and trimmed single-end fq
        #fastqc_raw_se = expand("fastqc/se/{sample}_fastqc.zip",
        #                        sample = SAMPLES["sample_id"]),
        mult_qc = "aggregated/QC_se/multiqc_report.html",
        return mult_qc + meth_qc + meta_quant + meth_filt


###############################
##  get corresponding bwa_index
def get_bwa_index():
    if config["spike_in"]:
        return REF.loc["bwa_idx_spikein"][1]
    else:
        return REF.loc["bwa_idx"][1]


################################################################
##  get raw fastq files for rename to consistant wildcard.sample
def get_raw_fastq(wildcards):
    if config["paired-end"]:
        R1 = SAMPLES.loc[wildcards.sample]["R1"],
        R2 = SAMPLES.loc[wildcards.sample]["R2"],
        return R1 + R2
    else:
        return SAMPLES.loc[wildcards.sample]["R1"]


#######################################################
##  get renamed fastq for FASQC and barcode extraction
def get_renamed_fastq(wildcards):
    if config["paired-end"]:
        R1 = "renamed_fq/{}_R1.fastq.gz".format(wildcards.sample),
        R2 = "renamed_fq/{}_R2.fastq.gz".format(wildcards.sample),
        return R1 + R2
    else:
        return "renamed_fq/{}.fastq.gz".format(wildcards.sample)


################################################
## get fastq for TRIM GALORE
## UMI extracted for paired-end reads, if exist
def get_fastq_4trim(wildcards):
    if config["paired-end"] == False:
        return "renamed_fq/{}.fastq.gz".format(wildcards.sample)
    elif config["add_umi"]:
        R1 = "barcoded_fq/{}_R1.fastq.gz".format(wildcards.sample),
        R2 = "barcoded_fq/{}_R2.fastq.gz".format(wildcards.sample),
        return R1 + R2
    else:
        R1 = "renamed_fq/{}_R1.fastq.gz".format(wildcards.sample),
        R2 = "renamed_fq/{}_R2.fastq.gz".format(wildcards.sample),
        return R1 + R2


##################################
## get trimmed fastq files for BWA
def get_trimmed_fastq(wildcards):
    if config["paired-end"]:
        R1 = "trimmed_fq/{}_R1_val_1.fq.gz".format(wildcards.sample),
        R2 = "trimmed_fq/{}_R2_val_2.fq.gz".format(wildcards.sample),
        return R1 + R2
    else:
        return "trimmed_fq/{}_trimmed.fq.gz".format(wildcards.sample)

########################################
## get dedup bam files for meth_qc_quant
def get_dedup_bam(wildcards):
    if config["paired-end"]:
        return "dedup_bam_pe/{}_dedup.bam".format(wildcards.sample)
    else:
        return "dedup_bam_se/{}_dedup.bam".format(wildcards.sample)


####################
## get FASTQC stats
def get_fastqc_stats():
    if config["paired-end"]:
        r1_raw  = expand("fastqc_pe/{samples}_R1_fastqc.zip", samples = SAMPLES["sample_id"]),
        r2_raw  = expand("fastqc_pe/{samples}_R2_fastqc.zip", samples = SAMPLES["sample_id"]),
        r1_trim = expand("fastqc_pe/{samples}_R1_val_1_fastqc.zip", samples = SAMPLES["sample_id"]),
        r2_trim = expand("fastqc_pe/{samples}_R2_val_2_fastqc.zip", samples = SAMPLES["sample_id"]),
        return r1_raw + r2_raw + r1_trim + r2_trim
    else:
         r1_raw  = expand("fastqc_se/{samples}_fastqc.zip", samples = SAMPLES["sample_id"]),
         r1_trim = expand("fastqc_se/{samples}_trimmed_fastqc.zip", samples = SAMPLES["sample_id"]),
         return r1_raw + r1_trim

######################
## get dedup bam stats
def get_dedup_bam_stats():
    if config["paired-end"]:
        return expand("dedup_bam_pe/{samples}_dedup.bam.stats.txt", samples = SAMPLES["sample_id"])
    else:
        return expand("dedup_bam_se/{samples}_dedup.bam.stats.txt", samples = SAMPLES["sample_id"])
