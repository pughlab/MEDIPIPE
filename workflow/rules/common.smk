## read in sample list
import pandas as pd

## read in sample and corresponding fq files talbe
SAMPLES = (
    pd.read_csv(config["samples"], sep="\t")
    .set_index("sample", drop=False)
    .sort_index()
)

## read in refrence files' info
REF = pd.read_csv(config["ref_files"], sep="\t", header = None, index_col = 0)

## working directory
wd = config["workdir"]
pipe_dir = config["pipeline_dir"]

#############################################
## get taget outputs based on the config file
#############################################

def get_rule_all_input():

    #map_out = expand("sorted_reads/{sample}_sorted.bam.bai", sample = SAMPLES["sample"]),
    dedup_out = expand("dedup_reads/{sample}_dedup.bam.bai", sample = SAMPLES["sample"]),
    medips_out = expand("medips/{sample}_qc_report.txt", sample = SAMPLES["sample"]),


    if config["paired-end"]:
        ## FASQC out for raw and trimmed paired-end fqs
        fastqc_raw_pe = expand("qc/pe/{sample}_{mate}_fastqc.html",
                                sample = SAMPLES["sample"],
                                mate = ["R1", "R2"]),
        fastqc_trimmed_r1 = expand("qc/pe/{sample}_R1_val_1_fastqc.html",
                                   sample = SAMPLES["sample"]),
        fastqc_trimmed_r2 = expand("qc/pe/{sample}_R2_val_2_fastqc.html",
                                   sample = SAMPLES["sample"]),
        fastqc_pe_out = fastqc_raw_pe + fastqc_trimmed_r1 + fastqc_trimmed_r2

        ## inferred insert size
        inferred_insert_size = expand("dedup_reads/{sample}_insert_size_histogram.pdf",
                                      sample = SAMPLES["sample"]),

        return fastqc_pe_out + dedup_out + inferred_insert_size + medips_out


    else:
        ## FASQC out for raw and trimmed single-end fq
        fastqc_raw_se = expand("qc/se/{sample}_fastqc.html",
                                sample = SAMPLES["sample"]),
        fastqc_trimmed_se = expand("qc/se/{sample}_trimmed_fastqc.html",
                                    sample = SAMPLES["sample"]),
        fastqc_se_out = fastqc_raw_se + fastqc_trimmed_se

        return fastqc_se_out + dedup_out + dedup_out + medips_out


###############################
##  get corresponding bwa_index
def get_bwa_index():
    if config["spike_in"]:
        return REF.loc["bwa_idx_spikein"][1]
    else:
        return REF.loc["bwa_idx"][1]


##############################################
##  get fastq files for FASTQC and TRIM_GALORE
def get_raw_fastq(wildcards):
    if config["paired-end"]:
        R1 = SAMPLES.loc[wildcards.sample]["R1"],
        R2 = SAMPLES.loc[wildcards.sample]["R2"],
        return R1 + R2
    else:
        return SAMPLES.loc[wildcards.sample]["R1"]


##################################
## get trimmed fastq files for BWA
def get_trimmed_fastq(wildcards):
    if config["paired-end"]:
        R1 = "trimmed_fq/{}_R1_val_1.fq.gz".format(wildcards.sample),
        R2 = "trimmed_fq/{}_R2_val_2.fq.gz".format(wildcards.sample),
        return R1 + R2
    else:
        return "trimmed_fq/{}_trimmed.fq.gz".format(wildcards.sample)
