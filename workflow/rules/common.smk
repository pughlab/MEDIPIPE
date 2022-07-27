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


## read in samples for aggregation
if config["aggreate"]:
    SAMPLES_AGGR = (
        pd.read_csv(config["samples_aggr"], sep="\t")
        .set_index("sample_id", drop=False)
        .sort_index()
    )
else:
    SAMPLES_AGGR = SAMPLES     ## must be defined

## paths for pipeline and/or reference data
wd = config["workdir"]
pipe_dir = config["pipeline_dir"]
umi_list = config["umi_list"]


## read in refrence files' info
REF = pd.read_csv(config["ref_files"], sep="\t", header = None, index_col = 0)
blacklist = REF.loc["blacklist"][1]   ## ENCODE blacklist


#############################################
## get taget outputs based on the config file
## either for individual samples or aggregate
## all samples listed in sample_aggr.tsv !!!
#############################################

def get_rule_all_input():
    ## ensure extra env installed
    extra_env = "extra_env/all_extra_env_installed",

    ## fixed outputs
    meth_qc = "aggregated/meth_qc.txt",
    meta_quant = "aggregated/meth_count.txt.gz",
    meth_filt = "autos_bfilt/meth_count_autos_bfilt.txt.gz",

    ######################################
    ## aggregated outputs for SAMPLES_aggr
    ## paired-end and spike-in
    if config["aggreate"] and config["paired-end"] and config["spike_in"]:
        mult_qc = "aggregated/QC_pe/multiqc_report.html",

        ## spike-ins
        spikein_mult_qc = "aggregated_spikein/QC_pe/multiqc_report.html",
        spikein_meth_qc = "aggregated_spikein/meth_qc.txt",
        spikein_meta_quant = "aggregated_spikein/meth_count.txt.gz",

        #fragment profiles
        fp_gc = "aggregated/fragment_profile_GC_corrected_1mb.tsv",        ## GC corrected fragment profile

        return  extra_env + mult_qc + meth_qc + meta_quant + meth_filt + spikein_mult_qc + spikein_meth_qc + spikein_meta_quant + fp_gc

    ## paired-end and no spike-in
    elif config["aggreate"] and config["paired-end"] and config["spike_in"] == False:
        mult_qc = "aggregated/QC_pe/multiqc_report.html",
        fp_gc = "aggregated/fragment_profile_GC_corrected_1mb.tsv",        ## GC corrected fragment profile

        return  extra_env + mult_qc + meth_qc + meta_quant + meth_filt + fp_gc

    ## single-end
    elif config["aggreate"] and config["paired-end"] == False:
        mult_qc = "aggregated/QC_se/multiqc_report.html",
        return extra_env + mult_qc + meth_qc + meta_quant + meth_filt


    #####################################
    ## outputs for each individual sample
    ## paired-end and spike-in
    elif config["aggreate"] == False and config["paired-end"] and config["spike_in"]:
        fq_qc = expand("fastqc_pe/{samples}_R2_fastqc.zip", samples = SAMPLES["sample_id"]),
        meth_out = expand("meth_qc_quant/{samples}_count.txt", samples = SAMPLES["sample_id"]),
        meth_spikein = expand("meth_qc_quant_spikein/{samples}_count.txt", samples = SAMPLES["sample_id"]),
        frag_size = expand("fragment_size/{samples}_insert_size_metrics.txt", samples = SAMPLES["sample_id"]),
        fp_gc = expand("fragment_profile/{samples}_50_Granges.bed", samples = SAMPLES["sample_id"]),

        return extra_env + fq_qc + frag_size + meth_out + meth_spikein + fp_gc

    ## paired-end without spike-in
    elif config["aggreate"] == False and config["paired-end"] and config["spike_in"] == False:
        fq_qc = expand("fastqc_pe/{samples}_R2_fastqc.zip", samples = SAMPLES["sample_id"]),
        meth_out = expand("meth_qc_quant/{samples}_count.txt", samples = SAMPLES["sample_id"]),
        frag_size = expand("fragment_size/{samples}_insert_size_metrics.txt", samples = SAMPLES["sample_id"]),
        fp_gc = expand("fragment_profile/{samples}_50_Granges.bed", samples = SAMPLES["sample_id"]),

        return extra_env + fq_qc + frag_size + meth_out + fp_gc

    ## single-end
    elif config["aggreate"] == False and config["paired-end"] == False:
         fq_qc = expand("fastqc_se/{samples}_R2_fastqc.zip", samples = SAMPLES["sample_id"]),
         meth_out = expand("meth_qc_quant/{samples}_count.txt", samples = SAMPLES["sample_id"]),

         return extra_env + fq_qc + meth_out




############################
## other functions for input
#############################

###############################
##  get corresponding bwa_index
def get_bwa_index():
    if config["spike_in"]:
        #return REF.loc["bwa_idx_spikein"][1]
        return config["spike_idx"]
    else:
        return REF.loc["bwa_index"][1]


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
    if config["paired-end"] and config["add_umi"]:
        R1 = "barcoded_fq_pe/{}_R1.fastq.gz".format(wildcards.sample),
        R2 = "barcoded_fq_pe/{}_R2.fastq.gz".format(wildcards.sample),
        return R1 + R2
    elif config["paired-end"] and config["add_umi"] == False:
        R1 = "renamed_fq/{}_R1.fastq.gz".format(wildcards.sample),
        R2 = "renamed_fq/{}_R2.fastq.gz".format(wildcards.sample),
        return R1 + R2
    elif config["paired-end"] == False and config["add_umi"]:
        return "barcoded_fq_se/{}.fastq.gz".format(wildcards.sample)
    else:
        return "renamed_fq/{}.fastq.gz".format(wildcards.sample)

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
    if config["paired-end"] and config["add_umi"]:
        return "dedup_bam_umi_pe/{}_dedup.bam".format(wildcards.sample)
    elif config["paired-end"] and config["add_umi"] == False:
        return "dedup_bam_pe/{}_dedup.bam".format(wildcards.sample)
    elif config["paired-end"] == False and config["add_umi"]:
        return "dedup_bam_umi_se/{}_dedup.bam".format(wildcards.sample)
    else:
        return "dedup_bam_se/{}_dedup.bam".format(wildcards.sample)


## spike-ins
def get_dedup_bam_spikein(wildcards):
    if config["paired-end"] and config["spike_in"]:
        return "dedup_bam_spikein/{}_spikein.bam".format(wildcards.sample)


#####################
## aggregaton #######
#####################

####################
## get FASTQC stats
def get_fastqc_stats():
    if config["paired-end"]:
        r1_raw  = expand("fastqc_pe/{samples}_R1_fastqc.zip", samples = SAMPLES_AGGR["sample_id"]),
        r2_raw  = expand("fastqc_pe/{samples}_R2_fastqc.zip", samples = SAMPLES_AGGR["sample_id"]),
        r1_trim = expand("fastqc_pe/{samples}_R1_val_1_fastqc.zip", samples = SAMPLES_AGGR["sample_id"]),
        r2_trim = expand("fastqc_pe/{samples}_R2_val_2_fastqc.zip", samples = SAMPLES_AGGR["sample_id"]),
        return r1_raw + r2_raw + r1_trim + r2_trim
    else:
         r1_raw  = expand("fastqc_se/{samples}_fastqc.zip", samples = SAMPLES_AGGR["sample_id"]),
         r1_trim = expand("fastqc_se/{samples}_trimmed_fastqc.zip", samples = SAMPLES_AGGR["sample_id"]),
         return r1_raw + r1_trim


######################
## get dedup bam stats
def get_dedup_bam_stats():
    if config["paired-end"] and config["add_umi"]:
        return expand("dedup_bam_umi_pe/{samples}_dedup.bam.stats.txt", samples = SAMPLES_AGGR["sample_id"])
    elif config["paired-end"] and config["add_umi"] == False:
        return expand("dedup_bam_pe/{samples}_dedup.bam.stats.txt", samples = SAMPLES_AGGR["sample_id"])
    elif config["paired-end"] == False and config["add_umi"]:
        return expand("dedup_bam_umi_se/{samples}_dedup.bam.stats.txt", samples = SAMPLES_AGGR["sample_id"])
    else:
        return expand("dedup_bam_se/{samples}_dedup.bam.stats.txt", samples = SAMPLES_AGGR["sample_id"])



#########################
## get spikeins bam stats
def get_spikein_stats():
    if config["spike_in"] and config["paired-end"]:
        bam_stats = expand("dedup_bam_spikein/{samples}_spikein.bam.stats.txt", samples = SAMPLES_AGGR["sample_id"]),
        frag_stats = expand("fragment_size_spikein/{samples}_insert_size_metrics.txt", samples = SAMPLES_AGGR["sample_id"]),
        return bam_stats + frag_stats
    else:
        return ""
