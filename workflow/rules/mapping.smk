## BWA alignment
rule bwa_map:
    input:
        #config["bwa_index"],
        get_bwa_index(),
        get_trimmed_fastq
    output:
        temp("mapped_reads/{sample}.bam")
    threads: 12
    log:
        "logs/{sample}_bwa_map.log"
    shell:
        "(bwa mem -M -t {threads}  {input} | "
        "samtools view -Sb --threads {threads} - > {output}) 2> {log}"

## fixmate, sort, index and stats bam file
rule samtools_sort_index_stats:
    input:
        "mapped_reads/{sample}.bam"
    output:
        bam = "sorted_reads/{sample}_sorted.bam",
        bai = "sorted_reads/{sample}_sorted.bam.bai",
        stat= "sorted_reads/{sample}_sorted.bam.stats.txt"
    threads: 12
    shell:
        ## --threads flag failed
        "(samtools fixmate -@ {threads} -m {input} - | "
        "samtools sort  -@ {threads} -o {output.bam} && "
        "samtools index -@ {threads} {output.bam} && "
        "samtools stats -@ {threads} {output.bam} > {output.stat})"

## markup, index and stats deduplicated file
rule samtools_markdup_stats:
    input:
        "sorted_reads/{sample}_sorted.bam"
    output:
        bam = "dedup_reads/{sample}_dedup.bam",
        bai = "dedup_reads/{sample}_dedup.bam.bai",
        stat= "dedup_reads/{sample}_dedup.bam.stats.txt"
    threads: 12
    shell:
        "(samtools markdup -@ {threads} -r {input} {output.bam} && "
        "samtools index -@ {threads} {output.bam} && "
        "samtools stats -@ {threads} {output.bam} > {output.stat})"

## infer insert size for paired-end reads_qc
rule insert_size:
    input:
        "dedup_reads/{sample}_dedup.bam"
    output:
        txt = "dedup_reads/{sample}_insert_size_metrics.txt",
        hist = "dedup_reads/{sample}_insert_size_histogram.pdf",
    log:
        "logs/{sample}_picard_insert_size.log"
    shell:
        "(java -jar /cluster/home/yzeng/miniconda3/envs/tcge-cfmedip-seq-pipeline/share/picard-2.26.6-0/picard.jar "
        "CollectInsertSizeMetrics M=0.05 I={input} O={output.txt} "
        "H={output.hist}) 2> {log}"
