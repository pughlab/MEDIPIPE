################################################################################
## multiQC summary
## using files instead of directories to ensure all samples qc metircs included
rule multiqc_pe:
    input:
        ## using "samples" to distinguish from wildcard.sample !!!
        get_fastqc_stats(),
        get_dedup_bam_stats(),
        expand("raw_bam/{samples}_sorted.bam.stats.txt",  samples = SAMPLES["sample_id"]),
        expand("dedup_bam_pe/{samples}_insert_size_metrics.txt", samples = SAMPLES["sample_id"]),
    output:
        # "summary/QC_se/multiqc_report.html"      ## only works for stand-alone mode,
        "summary/QC_pe/{sample}.html"              ## works for --cluster as we
    log:
        "logs/{sample}_pe.log"                        ## wildcard.sample needed for --cluster
    conda:
        "extra_env/multiQC.yaml"
    shell:
        "(multiqc {input} -o summary/QC_pe/) 2> {log}"


## single end
rule multiqc_se:
    input:
        ## using "samples" to distinguish from wildcard.sample !!!
        get_fastqc_stats(),
        get_dedup_bam_stats(),
        expand("raw_bam/{samples}_sorted.bam.stats.txt",  samples = SAMPLES["sample_id"]),
    output:
        "summary/QC_se/{sample}.html"
    log:
        "logs/{sample}_se.log"
    conda:
        "extra_env/multiQC.yaml"
    shell:
        "(multiqc {input} -o summary/QC_se/) 2> {log}"


##############################
## Aggregating meth QC reports
rule aggregate_meth_qc:
    input:
        ## using "samples" to distinguish from wildcard.sample !!!
        expand("meth_qc_quant/{samples}_meth_qc.txt", samples = SAMPLES["sample_id"])
    output:
        "summary/{sample}.txt"
    shell:
        "head -n 1 {input[0]} > {output} && "
        "cat {input} | sed '1~2d' >> {output}"


##########################################
## Aggregating meth Quantification outputs
## bin_id used to save storage space
rule aggregate_meth_quant:
    input:
        ## using "samples" to distinguish from wildcard.sample !!!
        bin  = expand("meth_qc_quant/{samples}_Granges_CpGs.bed", samples = SAMPLES["sample_id"][0]),
        cnt = expand("meth_qc_quant/{samples}_count.txt", samples = SAMPLES["sample_id"]),
        rpkm  = expand("meth_qc_quant/{samples}_rpkm.txt", samples = SAMPLES["sample_id"]),
        CNV_qsea   = expand("meth_qc_quant/{samples}_CNV_qsea.txt", samples = SAMPLES["sample_id"]),
        beta_qsea  = expand("meth_qc_quant/{samples}_beta_qsea.txt", samples = SAMPLES["sample_id"]),
        nrpm_qsea  = expand("meth_qc_quant/{samples}_nrpm_qsea.txt", samples = SAMPLES["sample_id"]),
        rms_medips = expand("meth_qc_quant/{samples}_rms_medips.txt", samples = SAMPLES["sample_id"]),
        rms_medestrand  = expand("meth_qc_quant/{samples}_rms_medestrand.txt", samples = SAMPLES["sample_id"]),
        logitbeta_qsea  = expand("meth_qc_quant/{samples}_logitbeta_qsea.txt", samples = SAMPLES["sample_id"]),
    output:
        bin  = "summary/{sample}_bin.bed",
        cnt  = "summary/{sample}_count.txt.gz",
        rpkm  = "summary/{sample}_rpkm.txt.gz",
        CNV_qsea   = "summary/{sample}_CNV_qseatxt.gz",
        beta_qsea  = "summary/{sample}_beta_qsea.txt.gz",
        nrpm_qsea  = "summary/{sample}_nrpm_qsea.txt.gz",
        rms_medips = "summary/{sample}_rms_medips.txt.gz",
        rms_medestrand  = "summary/{sample}_rms_medestrand.txt.gz",
        logitbeta_qsea  = "summary/{sample}_logitbeta_qsea.txt.gz"
    log:
        "logs/{sample}_quant_aggregate.log"
    resources:
        mem_mb=60000
    shell:
        "(cp {input.bin}  {output.bin} && "
        "paste {output.bin} {input.cnt}  | bgzip > {output.cnt} && tabix  {output.cnt} && "
        "paste {output.bin} {input.rpkm} | bgzip > {output.rpkm} && tabix  {output.rpkm} && "
        "paste {output.bin} {input.CNV_qsea}   |  bgzip > {output.CNV_qsea} && tabix  {output.CNV_qsea} && "
        "paste {output.bin} {input.beta_qsea}  |  bgzip > {output.beta_qsea} && tabix  {output.beta_qsea} && "
        "paste {output.bin} {input.nrpm_qsea}  |  bgzip > {output.nrpm_qsea} && tabix  {output.nrpm_qsea} && "
        "paste {output.bin} {input.rms_medips} |  bgzip > {output.rms_medips} && tabix  {output.rms_medips} && "
        "paste {output.bin} {input.rms_medestrand} | bgzip > {output.rms_medestrand} && tabix  {output.rms_medestrand} && "
        "paste {output.bin} {input.logitbeta_qsea} | bgzip > {output.logitbeta_qsea} && tabix  {output.logitbeta_qsea})  2> {log}"


##
"""
## filter out chrX, chrY, chrM and ENCODE blacklist bins
rule meth_bin_filter
    input:
        bin = "summary/{sample}_bin.bed"
    output:
        "summary/{sample}_autosomes_bin.bed",
        "summary/{sample}_autosomes_bfilt_bin.bed"
    params:
        blist = REF.loc["blacklist"][1]
    shell:
        "grep -v 'chrM\|chrX\|chrY' {input.bin}  >  {output[0]} && "
        "intersectBed -a {output} -b {params.blist} -v {output[1]} "

## meth quantification after filtering bins
rule meth_quant_filter
"""
