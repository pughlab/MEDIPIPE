## BWA alignment
rule bwa_map:
    input:
        config["bwa_index"],
        get_trimmed_fastq
    output:
        ("mapped_reads/{sample}.bam")
    threads: 2
    log:
        "logs/{sample}_bwa_map.log"
    shell:
        "(bwa mem -t {threads}  {input} | "
        "samtools view -Sb - > {output}) 2> {log}"

## fixmate, sort, index and stats bam file
rule samtools_sort_index:
    input:
        "mapped_reads/{sample}.bam"
    output:
        bam = "sorted_reads/{sample}_sorted.bam",
        bai = "sorted_reads/{sample}_sorted.bam.bai",
        stat= "sorted_reads/{sample}_sorted.bam.stats.txt"
    shell:
        "samtools fixmate -m {input} - | "
        "samtools sort -o {output.bam} && "
        "samtools index {output.bam} && "
        "samtools stats {output.bam} > {output.stat}"

## markup, index and stats deduplicated file
rule samtools_markdup:
    input:
        "sorted_reads/{sample}_sorted.bam"
    output:
        bam = "rmdup_reads/{sample}_rmdup.bam",
        bai = "rmdup_reads/{sample}_rmdup.bam.bai",
        stat= "rmdup_reads/{sample}_rmdup.bam.stats.txt"
    shell:
        "samtools markdup -r {input} {output.bam} && "
        "samtools index {output.bam} && "
        "samtools stats {output.bam} > {output.stat}"

## infer insert size for paired-end reads_qc
rule insert_size:
    input:
        "rmdup_reads/{sample}_rmdup.bam"
    output:
        txt = "rmdup_reads/{sample}_insert_size_metrics.txt",
        hist = "rmdup_reads/{sample}_insert_size_histogram.pdf",
    shell:
        "java -jar ~/miniconda3/envs/cfmedip-seq-pipeline/share/picard-2.26.6-0/picard.jar "
        "CollectInsertSizeMetrics M=0.05 I={input} O={output.txt} "
        "H={output.hist}"
