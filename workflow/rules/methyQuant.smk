## run MEDIPS for QC and quantifying
rule medips_medestrand:
    input:
        "dedup_bam/{sample}_dedup.bam"
    output:
        "methyQuant/{sample}_saturation.png",
        "methyQuant/{sample}_seqCoverage.png",
        "methyQuant/{sample}_qc_report.txt",
        "methyQuant/{sample}_meth_quant.txt",
        "methyQuant/{sample}_rpkm.wig",
        "methyQuant/{sample}_abs_methy.RDS"
    resources:
        mem_mb=60000
    params:
        ispaired = config["paired-end"],
        bsgenome = config["bsgenome"],
        scr_dir = config["pipeline_dir"],
    log:
        "logs/{sample}_medips_medestrand.log"
    conda:
        "conda_env_R.yaml"
    shell:
        "(Rscript --vanilla {params.scr_dir}/workflow/scripts/medips_medestrand.R "
        "{wildcards.sample} {input} {params.bsgenome} {params.ispaired} {output} "
        "{params.scr_dir}/workflow/dependencies/MeDEStrand) 2> {log}"
