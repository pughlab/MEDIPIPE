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
        bsgenome = config["bsgenome"]
    conda:                       ## failed to call env after created it in .snakemek/conda on iMAC
        "conda_env_R.yaml"
    shell:
        "Rscript --vanilla /cluster/home/yzeng/cfmedip-seq-pipeline/workflow/scripts/medips.R "
        "{wildcards.sample} {input} {params.bsgenome} {params.ispaired} {output}"


## run R script throuth script directive is still under testing
#script:
#    "/Users/yong/OneDrive_UHN/Projects/snakemake/cfmedip-seq-pipeline/workflow/rules/scripts/medips.R"
