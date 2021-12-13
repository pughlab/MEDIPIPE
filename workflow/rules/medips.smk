
## run MEDIPS for QC and quantifying
rule medips:
    input:
        #"rmdup_reads/{sample}_rmdup.bam"
        "rmdup_reads/test_in.txt"
    output:
        #"rmdup_reads/test_ori.txt",
        "rmdup_reads/test_out.txt",
    conda:                       ## failed to call env after created it in .snakemek/conda on iMAC
        "conda_env_R.yaml"
    shell:
        "Rscript --vanilla /Users/yong/OneDrive_UHN/Projects/snakemake/cfmedip-seq-pipeline/workflow/rules/scripts/medips.R {input} {output}"


## run R script throuth script directive is still under testing
#script:
#    "/Users/yong/OneDrive_UHN/Projects/snakemake/cfmedip-seq-pipeline/workflow/rules/scripts/medips.R"
