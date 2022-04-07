## Aggregating QC files
rule aggregate_qc:
    input:
        "raw_bam/{sample}_sorted.bam.stats.txt",
        "dedup_bam/{sample}_dedup.bam.stats.txt",
        "dedup_bam/{sample}_insert_size_metrics.txt",
        "meth_qc_quant/{sample}_meth_qc.txt"
    output:
        "aggregate_qc/bam_qc/{sample}_sorted.bam.stats.txt",
        "aggregate_qc/bam_qc/{sample}_dedup.bam.stats.txt",
        "aggregate_qc/bam_qc/{sample}_insert_size_metrics.txt",
        "aggregate_qc/meth_qc/{sample}_meth_qc.txt"
    shell:
        "cp {input[0]} {output[0]} && "
        "cp {input[1]} {output[1]} && "
        "cp {input[2]} {output[2]} && "
        "cp {input[3]} {output[3]}"


## Aggregating Quantification outputs
rule aggregate_quant:
    input:
        grange_cpg = expand("meth_qc_quant/{sample}_count.txt", sample = SAMPLES["sample_id"][1]),
        count = expand("meth_qc_quant/{sample}_count.txt", sample = SAMPLES["sample_id"]),
        rpkm  = expand("meth_qc_quant/{sample}_rpkm.txt", sample = SAMPLES["sample_id"]),
        rms_medips  = expand("meth_qc_quant/{sample}_rms_medips.txt", sample = SAMPLES["sample_id"]),
        rms_medestrand  = expand("meth_qc_quant/{sample}_rms_medestrand.txt", sample = SAMPLES["sample_id"]),
        CNV_qsea  = expand("meth_qc_quant/{sample}_CNV_qsea.txt", sample = SAMPLES["sample_id"]),
        beta_qsea  = expand("meth_qc_quant/{sample}_beta_qsea.txt", sample = SAMPLES["sample_id"]),
        nrpm_qsea  = expand("meth_qc_quant/{sample}_nrpm_qsea.txt", sample = SAMPLES["sample_id"]),
        logitbeta_qsea  = expand("meth_qc_quant/{sample}_logitbeta_qsea.txt", sample = SAMPLES["sample_id"]),
    output:
        count = "summary/count.txt.gz"
    shell:
        "cp {input.grange_cpg} {{output.grange_cpg} && "
        "paste {input.count} | gzip > {output.count}"
