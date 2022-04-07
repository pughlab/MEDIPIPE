## run MEDIPS + MeDEStrand + QSEA for QC and quantification
## relative/absolute methylation levels will be estimated

rule meth_qc_quant:
    input:
        "dedup_bam/{sample}_dedup.bam"
    output:
        "meth_qc_quant/{sample}_saturation.png",
        "meth_qc_quant/{sample}_seqCoverage.png",
        "meth_qc_quant/{sample}_meth_qc.txt",
        "meth_qc_quant/{sample}_meth_quant.RData"
    resources:
        mem_mb=60000
    params:
        ispaired = config["paired-end"],
        bsgenome = config["bsgenome"],
        scr_dir = config["pipeline_dir"],
    log:
        "logs/{sample}_medips_medestrand_qsea.log"
    conda:
        "extra_env/conda_env_R.yaml"
    shell:
        "(Rscript --vanilla {params.scr_dir}/workflow/scripts/medips_medestrand_qsea.R "
        "{wildcards.sample} {input} {params.bsgenome} {params.ispaired} "
        "{params.scr_dir}/workflow/dependencies/MeDEStrand) 2> {log}"
