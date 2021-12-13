





     clean(r'''
        samtools view -buS -f 2 -F 4 -@32 {input} |
        samtools fixmate -m - - |
        samtools sort -@32 -o {output.bam} && samtools index {output.bam}
        ''')


    "samtools sort -T sorted_reads/{wildcards.sample} "
    "-O bam {input} > {output}"

    ## add index
    rule samtools_index:
        input:
            "sorted_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam.bai"
        shell:
            "samtools index {input}"



~/OneDrive_UHN/Projects/snakemake/cfmedip-seq-pipeline

snakemake -p --snakefile ./workflow/Snakefile --configfile ./workflow/config/config.yaml

/Users/yong/tmp/sorted_reads
