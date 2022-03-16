## Aggregating QC files
rule aggregate_qc:
    input:
        "raw_bam/{sample}_sorted.bam.stats.txt"
        "dedup_bam/{sample}_dedup.bam.stats.txt"
        "dedup_bam/{sample}_insert_size_metrics.txt",
        "meth_quant/{sample}_meth_qc_report.txt",
    output:
        "aggregate_qc/bam_qc/{sample}_sorted.bam.stats.txt"
        "aggregate_qc/bam_qc/{sample}_dedup.bam.stats.txt"
        "aggregate_qc/bam_qc/{sample}_insert_size_metrics.txt",
        "aggregate_qc/meth_qc/{sample}_meth_qc_report.txt",

    shell:
        "cp {input[0]} {output[0]} && "
        "cp {input[1]} {output[1]} && "
        "cp {input[2]} {output[2]} && "
        "cp {input[3]} {output[3]}"
