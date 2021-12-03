rule bwa_map:
    input:
        config["bwa_index"],
        "tmp/trimmed/{sample}_trimmed.fq"
    output:
        temp("tmp/mapped_reads/{sample}.bam")
    threads: 2
    log:
        "tmp/logs/{sample}_bwa_map.log"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    shell:
        "(bwa mem -R '{params.rg}' -t {threads}  {input} | "
        "samtools view -Sb - > {output}) 2> {log}"


rule samtools_sort:
    input:
        "tmp/mapped_reads/{sample}.bam"
    output:
        "tmp/sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T tmp/sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        "tmp/sorted_reads/{sample}.bam"
    output:
        "tmp/sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"
