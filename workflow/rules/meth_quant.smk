## run MEDIPS for QC and quantifying
rule medips_medestrand:
    input:
        "dedup_bam/{sample}_dedup.bam"
    output:
        "meth_quant/{sample}_saturation.png",
        "meth_quant/{sample}_seqCoverage.png",
        "meth_quant/{sample}_meth_qc_report.txt",
        "meth_quant/{sample}_meth_quant.txt",
        "meth_quant/{sample}_meth_rpkm.wig",
        "meth_quant/{sample}_meth_abs.RDS"
    resources:
        mem_mb=60000
    params:
        ispaired = config["paired-end"],
        bsgenome = config["bsgenome"],
        scr_dir = config["pipeline_dir"],
    log:
        "logs/{sample}_medips_medestrand.log"
    conda:
        "extra_env/conda_env_R.yaml"
    shell:
        "(Rscript --vanilla {params.scr_dir}/workflow/scripts/medips_medestrand.R "
        "{wildcards.sample} {input} {params.bsgenome} {params.ispaired} {output} "
        "{params.scr_dir}/workflow/dependencies/MeDEStrand) 2> {log}"
