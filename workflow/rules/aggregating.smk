################################################################################
## multiQC
## using files instead of directories to ensure all samples qc metircs included
## require sample_aggr.tsv !!!
################################################################################
rule multiqc_pe:
    input:
        ## using "samples" to distinguish from wildcard.sample !!!
        get_fastqc_stats(),
        get_dedup_bam_stats(),
        expand("raw_bam/{samples}_sorted.bam.stats.txt",  samples = SAMPLES_AGGR["sample_id"]),
        expand("fragment_size/{samples}_insert_size_metrics.txt", samples = SAMPLES_AGGR["sample_id"]),
    output:
        # "aggregated/QC_se/multiqc_report.html"      ## only works for stand-alone mode,
        "aggregated/QC_pe/{sample}_report.html",
        ## inputs for aggr_qc_report
        "aggregated/QC_pe/multiqc_data/{sample}_fastqc.txt",
        "aggregated/QC_pe/multiqc_data/{sample}_samtools_stats.txt",
        "aggregated/QC_pe/multiqc_data/{sample}_picard_insertSize.txt"
    log:
        "logs/{sample}_pe.log"                        ## wildcard.sample needed for --cluster
    conda:
        "extra_env/multiQC.yaml"
    shell:
        "(multiqc -f {input} -o aggregated/QC_pe/) 2> {log}"


## QC for sipke-ins separately
rule multiqc_pe_spikein:
    input:
        get_spikein_stats()
    output:
        # "aggregated/QC_se/multiqc_report.html"      ## only works for stand-alone mode,
        "aggregated_spikein/QC_pe/{sample}.html"              ## works for --cluster as we
    log:
        "logs/{sample}_pe_spikein.log"                        ## wildcard.sample needed for --cluster
    conda:
        "extra_env/multiQC.yaml"
    shell:
        "(multiqc -f {input} -o aggregated_spikein/QC_pe/) 2> {log}"


## single end
rule multiqc_se:
    input:
        ## using "samples" to distinguish from wildcard.sample !!!
        get_fastqc_stats(),
        get_dedup_bam_stats(),
        expand("raw_bam/{samples}_sorted.bam.stats.txt",  samples = SAMPLES_AGGR["sample_id"]),
    output:
        "aggregated/QC_se/{sample}.html",
        ## inputs for aggr_qc_report
        "aggregated/QC_se/multiqc_data/{sample}_fastqc.txt",
        "aggregated/QC_se/multiqc_data/{sample}_samtools_stats.txt"
    log:
        "logs/{sample}_se.log"
    conda:
        "extra_env/multiQC.yaml"
    shell:
        "(multiqc -f {input} -o aggregated/QC_se/) 2> {log}"


##############################
## Aggregating meth QC reports
##############################
rule aggregate_meth_qc:
    input:
        ## using "samples" to distinguish from wildcard.sample !!!
        expand("meth_qc_quant/{samples}_meth_qc.txt", samples = SAMPLES_AGGR["sample_id"])
    output:
        "aggregated/{sample}_qc.txt"
    shell:
        "head -n 1 {input[0]} > {output} && "
        "cat {input} | sed '1~2d' >> {output}"

## Aggregating spike-ins QC
rule aggregate_meth_qc_spikein:
    input:
        ## using "samples" to distinguish from wildcard.sample !!!
        expand("meth_qc_quant_spikein/{samples}_meth_qc.txt", samples = SAMPLES_AGGR["sample_id"])
    output:
        "aggregated_spikein/{sample}.txt"
    shell:
        "head -n 1 {input[0]} > {output} && "
        "cat {input} | sed '1~2d' >> {output}"


###############################################
## Aggregating QC reports: matrix + html report
###############################################
rule aggregate_qc_report:
    input:
        get_aggr_qc_stats()
    output:
        "aggregated/{sample}_qc_report.html"
    resources:
        mem_mb=60000
    params:
        scr_dir = config["pipe_dir"],
        out_dir = config['work_dir'],
        aggr = config["samples_aggr"],
        ispaired = config["paired-end"],
    log:
        "logs/{sample}_qc_report.log"
    conda:
        "extra_env/R_aggr_qc_report.yaml"
    shell:
        "(Rscript --vanilla {params.scr_dir}/workflow/scripts/aggr_qc_report.R "
        "{params.aggr} {params.ispaired} {input} {params.scr_dir} {params.out_dir}) 2> {log}"


################################
## Aggregating fragment profiles
################################
rule aggregate_fragment_profile:
    input:
        bin1 = expand("fragment_profile/{samples}_10_Granges.bed", samples = SAMPLES_AGGR["sample_id"][0]),
        bin5 = expand("fragment_profile/{samples}_50_Granges.bed", samples = SAMPLES_AGGR["sample_id"][0]),
        mb1 = expand("fragment_profile/{samples}_10_100kb_fragment_profile_GC_corrected_Ratio.txt", samples = SAMPLES_AGGR["sample_id"]),
        mb5 = expand("fragment_profile/{samples}_50_100kb_fragment_profile_GC_corrected_Ratio.txt", samples = SAMPLES_AGGR["sample_id"])
    output:
        bin1 = "aggregated/{sample}_Granges_1mb.bed",
        bin5 = "aggregated/{sample}_Granges_5mb.bed",
        mb1 = "aggregated/{sample}_GC_corrected_1mb.tsv",
        mb5 = "aggregated/{sample}_GC_corrected_5mb.tsv",
    log:
        "logs/{sample}_aggregate.log"
    resources:
        mem_mb=60000
    shell:
        "(cp {input.bin1} {output.bin1} && paste {output.bin1} {input.mb1}  >  {output.mb1} && "
        " cp {input.bin5} {output.bin5} && paste {output.bin5} {input.mb5}  >  {output.mb5}) 2> {log}"


##########################################
## Aggregating meth Quantification outputs
## bin_id used to save storage space
##########################################
rule aggregate_meth_quant:
    input:
        ## using "samples" to distinguish from wildcard.sample !!!
        bin  = expand("meth_qc_quant/{samples}_Granges_CpGs.bed", samples = SAMPLES_AGGR["sample_id"][0]),
        cnt = expand("meth_qc_quant/{samples}_count.txt", samples = SAMPLES_AGGR["sample_id"]),
        rpkm  = expand("meth_qc_quant/{samples}_rpkm.txt", samples = SAMPLES_AGGR["sample_id"]),
        CNV_qsea   = expand("meth_qc_quant/{samples}_CNV_qsea.txt", samples = SAMPLES_AGGR["sample_id"]),
        beta_qsea  = expand("meth_qc_quant/{samples}_beta_qsea.txt", samples = SAMPLES_AGGR["sample_id"]),
        nrpm_qsea  = expand("meth_qc_quant/{samples}_nrpm_qsea.txt", samples = SAMPLES_AGGR["sample_id"]),
        rms_medips = expand("meth_qc_quant/{samples}_rms_medips.txt", samples = SAMPLES_AGGR["sample_id"]),
        rms_medestrand  = expand("meth_qc_quant/{samples}_rms_medestrand.txt", samples = SAMPLES_AGGR["sample_id"]),
        logitbeta_qsea  = expand("meth_qc_quant/{samples}_logitbeta_qsea.txt", samples = SAMPLES_AGGR["sample_id"]),
    output:
        bin  = "aggregated/{sample}_bin.bed",
        cnt  = "aggregated/{sample}_count.txt.gz",
        rpkm  = "aggregated/{sample}_rpkm.txt.gz",
        CNV_qsea   = "aggregated/{sample}_CNV_qsea.txt.gz",
        beta_qsea  = "aggregated/{sample}_beta_qsea.txt.gz",
        nrpm_qsea  = "aggregated/{sample}_nrpm_qsea.txt.gz",
        rms_medips = "aggregated/{sample}_rms_medips.txt.gz",
        rms_medestrand  = "aggregated/{sample}_rms_medestrand.txt.gz",
        logitbeta_qsea  = "aggregated/{sample}_logitbeta_qsea.txt.gz"
    log:
        "logs/{sample}_quant_aggregate.log"
    resources:
        mem_mb=60000
    shell:
        "(cp {input.bin}  {output.bin} && "
        "paste {output.bin} {input.cnt}  | bgzip > {output.cnt} && tabix -p bed {output.cnt} && "
        "paste {output.bin} {input.rpkm} | bgzip > {output.rpkm} && tabix -p bed {output.rpkm} && "
        "paste {output.bin} {input.CNV_qsea}   |  bgzip > {output.CNV_qsea} && tabix -p bed {output.CNV_qsea} && "
        "paste {output.bin} {input.beta_qsea}  |  bgzip > {output.beta_qsea} && tabix -p bed {output.beta_qsea} && "
        "paste {output.bin} {input.nrpm_qsea}  |  bgzip > {output.nrpm_qsea} && tabix -p bed {output.nrpm_qsea} && "
        "paste {output.bin} {input.rms_medips} |  bgzip > {output.rms_medips} && tabix -p bed {output.rms_medips} && "
        "paste {output.bin} {input.rms_medestrand} | bgzip > {output.rms_medestrand} && tabix -p bed {output.rms_medestrand} && "
        "paste {output.bin} {input.logitbeta_qsea} | bgzip > {output.logitbeta_qsea} && tabix -p bed {output.logitbeta_qsea})  2> {log}"


## aggregate spike-ins
rule aggregate_meth_quant_spikein:
    input:
        ## using "samples" to distinguish from wildcard.sample !!!
        bin  = expand("meth_qc_quant_spikein/{samples}_Granges_CpGs.bed", samples = SAMPLES_AGGR["sample_id"][0]),
        cnt = expand("meth_qc_quant_spikein/{samples}_count.txt", samples = SAMPLES_AGGR["sample_id"]),
        rpkm  = expand("meth_qc_quant_spikein/{samples}_rpkm.txt", samples = SAMPLES_AGGR["sample_id"]),
    output:
        bin  = "aggregated_spikein/{sample}_bin.bed",
        cnt  = "aggregated_spikein/{sample}_count.txt.gz",
        rpkm  = "aggregated_spikein/{sample}_rpkm.txt.gz",
    log:
        "logs/{sample}_quant_aggregate_spikein.log"
    resources:
        mem_mb=60000
    shell:
        "(cp {input.bin}  {output.bin} && "
        "paste {output.bin} {input.cnt}  | bgzip > {output.cnt} && tabix -p bed {output.cnt} && "
        "paste {output.bin} {input.rpkm} | bgzip > {output.rpkm} && tabix -p bed {output.rpkm})  2> {log}"



########################################################
## filter out chrX, chrY, chrM and ENCODE blacklist bins
## autos_bfilt: autosomes + blacklist fitered
########################################################
rule meth_bin_filter:
    input:
        bin = "aggregated/{sample}_bin.bed"
    output:
        "autos_bfilt/{sample}_autos_bfilt_bin.bed",
        "autos_bfilt/{sample}_autos_bfilt_bin_merged.bed"
    conda:
        "extra_env/bedtools.yaml"
    params:
        blacklist
    shell:
        ## autosomes
        "head -1 {input.bin} > autos_bfilt/tmp_header.bed && "
        "grep -v 'chrM\|chrX\|chrY' {input.bin}  >  autos_bfilt/tmp_1.bed && "
        ## mask blacklist
        "intersectBed -a autos_bfilt/tmp_1.bed -b {params} -v >  autos_bfilt/tmp_2.bed && "
        "sort -k4,4n autos_bfilt/tmp_2.bed > autos_bfilt/tmp_3.bed && "
        "cat autos_bfilt/tmp_header.bed  autos_bfilt/tmp_3.bed > {output[0]} && "

        ## merge filtered bins for tabix
        "bedtools merge -i {output[0]} -d 1 | sort -V -k1,1 -k2,2n > {output[1]} && "
        "rm autos_bfilt/tmp_*.bed "


#################################################################################
## meth quantification after filtering chrX, chrY, chrM and ENCODE blacklist bins
rule meth_quant_filter:
    input:
        bin = "autos_bfilt/{sample}_autos_bfilt_bin_merged.bed",
        cnt = "aggregated/{sample}_count.txt.gz",
        rpkm  = "aggregated/{sample}_rpkm.txt.gz",
        CNV_qsea   = "aggregated/{sample}_CNV_qsea.txt.gz",
        beta_qsea  = "aggregated/{sample}_beta_qsea.txt.gz",
        nrpm_qsea  = "aggregated/{sample}_nrpm_qsea.txt.gz",
        rms_medips = "aggregated/{sample}_rms_medips.txt.gz",
        rms_medestrand  = "aggregated/{sample}_rms_medestrand.txt.gz",
        logitbeta_qsea  = "aggregated/{sample}_logitbeta_qsea.txt.gz",
    output:
        cnt = "autos_bfilt/{sample}_count_autos_bfilt.txt.gz",
        rpkm  = "autos_bfilt/{sample}_rpkm_autos_bfilt.txt.gz",
        CNV_qsea   = "autos_bfilt/{sample}_CNV_qsea_autos_bfilt.txt.gz",
        beta_qsea  = "autos_bfilt/{sample}_beta_qsea_autos_bfilt.txt.gz",
        nrpm_qsea  = "autos_bfilt/{sample}_nrpm_qsea_autos_bfilt.txt.gz",
        rms_medips = "autos_bfilt/{sample}_rms_medips_autos_bfilt.txt.gz",
        rms_medestrand  = "autos_bfilt/{sample}_rms_medestrand_autos_bfilt.txt.gz",
        logitbeta_qsea  = "autos_bfilt/{sample}_logitbeta_qsea_autos_bfilt.txt.gz"
    resources:
        mem_mb=60000
    shell:
        "tabix -p bed -R {input.bin} -h {input.cnt} | bgzip > {output.cnt} && tabix -p bed {output.cnt} && "
        "tabix -p bed -R {input.bin} -h {input.rpkm} | bgzip > {output.rpkm} && tabix -p bed {output.rpkm} && "
        "tabix -p bed -R {input.bin} -h {input.CNV_qsea}   |  bgzip > {output.CNV_qsea} && tabix -p bed {output.CNV_qsea} && "
        "tabix -p bed -R {input.bin} -h {input.beta_qsea}  |  bgzip > {output.beta_qsea} && tabix -p bed {output.beta_qsea} && "
        "tabix -p bed -R {input.bin} -h {input.nrpm_qsea}  |  bgzip > {output.nrpm_qsea} && tabix -p bed {output.nrpm_qsea} && "
        "tabix -p bed -R {input.bin} -h {input.rms_medips} |  bgzip > {output.rms_medips} && tabix -p bed {output.rms_medips} && "
        "tabix -p bed -R {input.bin} -h {input.rms_medestrand} | bgzip > {output.rms_medestrand} && tabix -p bed {output.rms_medestrand} && "
        "tabix -p bed -R {input.bin} -h {input.logitbeta_qsea} | bgzip > {output.logitbeta_qsea} && tabix -p bed {output.logitbeta_qsea}"
