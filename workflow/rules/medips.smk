## run MEDIPS for QC and quantifying
rule medips:
    input:
        "dedup_reads/{sample}_dedup.bam"
    output:
        "medips/{sample}_saturation.png",
        "medips/{sample}_seqCoverage.png",
        "medips/{sample}_qc_report.txt",
        "medips/{sample}_meth_quant.txt",
        "medips/{sample}_rpkm.wig"
    params:
        ispaired = config["paired-end"],
        bsgenome = config["bsgenome"],
        scr_dir = config["pipeline_dir"]
    log:
        "logs/{sample}_medips.log"
    conda:
        "conda_env_R.yaml"
    shell:
        "(Rscript --vanilla {params.scr_dir}/workflow/scripts/medips.R "
        "{wildcards.sample} {input} {params.bsgenome} {params.ispaired} {output}) 2> {log}"


## run R script throuth script directive is still under testing
#script:
#    "/Users/yong/OneDrive_UHN/Projects/snakemake/cfmedip-seq-pipeline/workflow/rules/scripts/medips.R"
