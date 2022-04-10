###########################################################
## run MEDIPS + MeDEStrand + QSEA for QC and quantification
## relative/absolute methylation levels will be estimated
rule meth_qc_quant:
    input:
        get_dedup_bam
    output:
        "meth_qc_quant/{sample}_meth_qc.txt",
        "meth_qc_quant/{sample}_meth_quant.RData",
        temp("meth_qc_quant/{sample}_Granges_CpGs.bed"),
        "meth_qc_quant/{sample}_count.txt",
        "meth_qc_quant/{sample}_rpkm.txt",
        "meth_qc_quant/{sample}_CNV_qsea.txt",
        "meth_qc_quant/{sample}_beta_qsea.txt",
        "meth_qc_quant/{sample}_nrpm_qsea.txt",
        "meth_qc_quant/{sample}_rms_medips.txt",
        "meth_qc_quant/{sample}_rms_medestrand.txt",
        "meth_qc_quant/{sample}_logitbeta_qsea.txt"
    resources:
        mem_mb=60000
    params:
        ispaired = config["paired-end"],
        bsgenome = config["bsgenome"],
        scr_dir = config["pipeline_dir"],
    log:
        "logs/{sample}_medips_medestrand_qsea.log"
    conda:
        "extra_env/R.yaml"
    shell:
        "(Rscript --vanilla {params.scr_dir}/workflow/scripts/medips_medestrand_qsea.R "
        "{wildcards.sample} {input} {params.bsgenome} {params.ispaired} "
        "{params.scr_dir}/workflow/dependencies/MeDEStrand) 2> {log}"
