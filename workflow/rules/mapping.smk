################
## BWA alignment
rule bwa_map:
    input:
        #config["bwa_index"],
        get_bwa_index(),
        get_trimmed_fastq
    output:
        temp("raw_bam/{sample}.bam")
    threads: 12
    log:
        "logs/{sample}_bwa_map.log"
    shell:
        "(bwa mem -M -t {threads}  {input} | "
        "samtools view -Sb --threads {threads} - > {output}) 2> {log}"


##########################################
## raw bams without any filtering
## fixmate, sort, index and stats bam file
rule samtools_sort_index_stats:
    input:
        "raw_bam/{sample}.bam"
    output:
        bam = "raw_bam/{sample}_sorted.bam",
        #bai = "raw_bam/{sample}_sorted.bam.bai",
        stat= "raw_bam/{sample}_sorted.bam.stats.txt"
    threads: 12
    shell:
        ## --threads flag failed
        "(samtools fixmate -@ {threads} -m {input} - | "
        "samtools sort  -@ {threads} -o {output.bam} && "
        "samtools index -@ {threads} {output.bam} && "
        "samtools stats -@ {threads} {output.bam} > {output.stat})"


##########################################################################
## to filter out unmapped & non-uniquely mapped, not properly paired reads
## Deduplication with markup, index and stats deduplicated file
## No UMIs !!!!
rule samtools_markdup_stats:
    input:
        "raw_bam/{sample}_sorted.bam"
    output:
        bam = "dedup_bam_pe/{sample}_dedup.bam",
        #bai = "dedup_bam/{sample}_dedup.bam.bai",
        stat= "dedup_bam_pe/{sample}_dedup.bam.stats.txt"
    threads: 12
    shell:
        "(samtools view -b -f 2 -F 2828 --threads {threads} {input} | "
        "samtools markdup -@ {threads} -r - {output.bam} && "
        "samtools index -@ {threads} {output.bam} && "
        "samtools stats -@ {threads} {output.bam} > {output.stat})"

"""
## to filter out unmapped & non-uniquely mapped, not properly paired reads
## Deduplication with UMI-tools, which takes UMI and coordinates info into account
## UMI-tools dosen't support parallel threads yet!!
## index and stats deduplicated bams
rule samtools_umi_tools:
    input:
        "raw_bam/{sample}_sorted.bam"
    output:
        temp(filter_bam = "dedup_bam_pe_umi/{sample}_filter.bam"),
        dedup_bam = "dedup_bam_pe_umi/{sample}_dedup.bam",
        umi_stat = "dedup_bam_pe_umi/{sample}_UMI.stats.txt",
        bam_stat= "dedup_bam_pe_umi/{sample}_dedup.bam.stats.txt"
    log:
        "logs/{sample}_dedup_umi.log"
    conda:
        "extra_env/UMI-tools.yaml"
    threads: 12
    shell:
        "(samtools view -b -f 2 -F 2828 --threads {threads} {input} > {output.filter_bam} && "
        "umi_tools dedup --paired -I {output.filter_bam} -S {output.dedup_bam} --output-stats={output.umi_stat} && "
        "samtools index -@ {threads} {output.dedup_bam} && "
        "samtools stats -@ {threads} {output.dedup_bam} > {output.bam_stat}) 2> {log}"

"""

## extract spike-ins bam after deduplication
## paired-end only so far
rule samtools_spikein_sort_index_stats:
    input:
        #"raw_bam/{sample}_sorted.bam"      ## lead to ambiguous wildcards!?
        "dedup_bam_pe/{sample}_dedup.bam"
    output:
        bam = "dedup_bam_spikein/{sample}_spikein.bam",
        #bai = "raw_bam/{sample}_sorted.bam.bai",
        stat= "dedup_bam_spikein/{sample}_spikein.bam.stats.txt"
    threads: 12
    params:
        spikein_chr = config["spike_in_chr"]
    shell:
        ## --threads flag failed
        "(samtools view  -@ {threads} -hbS {input} {params.spikein_chr} | "
        "samtools  sort  -@ {threads} -o {output.bam} && "
        "samtools  index -@ {threads} {output.bam} && "
        "samtools  stats -@ {threads} {output.bam} > {output.stat})"


## single-end
## filtering differ
rule samtools_markdup_stats_se:
    input:
        "raw_bam/{sample}_sorted.bam"
    output:
        bam = "dedup_bam_se/{sample}_dedup.bam",
        #bai = "dedup_bam/{sample}_dedup.bam.bai",
        stat= "dedup_bam_se/{sample}_dedup.bam.stats.txt"
    threads: 12
    shell:
        "(samtools view -b -F 2820 --threads {threads} {input} | "
        "samtools markdup -@ {threads} -r - {output.bam} && "
        "samtools index -@ {threads} {output.bam} && "
        "samtools stats -@ {threads} {output.bam} > {output.stat})"

############################################
## infer insert size for paired-end reads_qc
rule insert_size:
    input:
        "dedup_bam_pe/{sample}_dedup.bam"
    output:
        txt = "dedup_bam_pe/{sample}_insert_size_metrics.txt",
        hist = "dedup_bam_pe/{sample}_insert_size_histogram.pdf"
    params:
        pipeline_env = config["pipeline_env"]
    log:
        "logs/{sample}_picard_insert_size.log"
    shell:
        "(java -jar {params.pipeline_env}/share/picard-2.26.6-0/picard.jar "
        "CollectInsertSizeMetrics M=0.05 I={input} O={output.txt} "
        "H={output.hist}) 2> {log}"

############################################
## infer insert size for paired-end reads_qc
## spike-ins
rule insert_size_spikein:
    input:
        "dedup_bam_spikein/{sample}_spikein.bam"
    output:
        txt = "dedup_bam_spikein/{sample}_insert_size_metrics.txt",
        hist = "dedup_bam_spikein/{sample}_insert_size_histogram.pdf"
    params:
        pipeline_env = config["pipeline_env"]
    log:
        "logs/{sample}_picard_insert_size_spikein.log"
    shell:
        "(java -jar {params.pipeline_env}/share/picard-2.26.6-0/picard.jar "
        "CollectInsertSizeMetrics M=0.05 I={input} O={output.txt} "
        "H={output.hist}) 2> {log}"
