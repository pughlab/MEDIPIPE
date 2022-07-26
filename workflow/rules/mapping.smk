################
## BWA alignment
################
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
###########################################################################

###############
## without UMIs
rule samtools_markdup_stats_pe:
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

## single-end filtering differ
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


#######################################################################################
## with UMIs !!!!
## Deduplication with UMI-tools, which takes both UMI and coordinates info into account
## UMI-tools dosen't support parallel threads yet!!
## paired-end
rule samtools_umi_tools_pe:
    input:
        "raw_bam/{sample}_sorted.bam"
    output:
        dedup_bam = "dedup_bam_umi_pe/{sample}_dedup.bam",
        bam_stat = "dedup_bam_umi_pe/{sample}_dedup.bam.stats.txt",
    params:
        tmp_bam = "dedup_bam_umi_pe/{sample}_tmp.bam",
        stat_prefix = "dedup_bam_umi_pe/{sample}_dedup"
    threads: 12
    conda:
        "extra_env/umi_tools.yaml"
    log:
        "logs/{sample}_dedup_umi.log"
    shell:
        "(umi_tools dedup --paired -I {input} -S {params.tmp_bam} --umi-separator=':' --output-stats={params.stat_prefix} && "
        "samtools view -b -F 2820 --threads {threads} {params.tmp_bam} > {output.dedup_bam} && "
        "samtools index -@ {threads} {output.dedup_bam}  && rm {params.tmp_bam} && "
        "samtools stats -@ {threads} {output.dedup_bam} > {output.bam_stat}) 2> {log}"


## single-end with UMIs: different samtools filtering flags
rule samtools_umi_tools_se:
    input:
        "raw_bam/{sample}_sorted.bam"
    output:
        dedup_bam = "dedup_bam_umi_se/{sample}_dedup.bam",
        bam_stat = "dedup_bam_umi_se/{sample}_dedup.bam.stats.txt",
    params:
        tmp_bam = "dedup_bam_umi_se/{sample}_tmp.bam",
        stat_prefix = "dedup_bam_umi_se/{sample}_dedup"
    threads: 12
    conda:
        "extra_env/umi_tools.yaml"
    log:
        "logs/{sample}_dedup_umi.log"
    shell:
        "(umi_tools dedup --paired -I {input} -S {params.tmp_bam} --umi-separator=':' --output-stats={params.stat_prefix} && "
        "samtools view -b -f 2 -F 2828 --threads {threads} {params.tmp_bam} > {output.dedup_bam} && "
        "samtools index -@ {threads} {output.dedup_bam}  && rm {params.tmp_bam} && "
        "samtools stats -@ {threads} {output.dedup_bam} > {output.bam_stat}) 2> {log}"



############################################
## extract spike-ins bam after deduplication
############################################
## paired-end only so far !!
rule samtools_spikein_sort_index_stats:
    input:
        #"raw_bam/{sample}_sorted.bam"      ## lead to ambiguous wildcards!?
        #"dedup_bam_pe/{sample}_dedup.bam"
        get_dedup_bam
    output:
        bam = "dedup_bam_spikein/{sample}_spikein.bam",
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



############################################
## infer insert size for paired-end reads_qc
############################################

rule insert_size:
    input:
        #"dedup_bam_pe/{sample}_dedup.bam"
        get_dedup_bam
    output:
        txt = "fragment_size/{sample}_insert_size_metrics.txt",
        hist = "fragment_size/{sample}_insert_size_histogram.pdf"
    params:
        pipeline_env = config["pipeline_env"]
    log:
        "logs/{sample}_picard_insert_size.log"
    shell:
        "(java -jar {params.pipeline_env}/share/picard-2.26.6-0/picard.jar "
        "CollectInsertSizeMetrics M=0.05 I={input} O={output.txt} "
        "H={output.hist}) 2> {log}"


## spike-ins

rule insert_size_spikein:
    input:
        "dedup_bam_spikein/{sample}_spikein.bam"
    output:
        txt = "fragment_size_spikein/{sample}_insert_size_metrics.txt",
        hist = "fragment_size_spikein/{sample}_insert_size_histogram.pdf"
    params:
        pipeline_env = config["pipeline_env"]
    log:
        "logs/{sample}_picard_insert_size_spikein.log"
    shell:
        "(java -jar {params.pipeline_env}/share/picard-2.26.6-0/picard.jar "
        "CollectInsertSizeMetrics M=0.05 I={input} O={output.txt} "
        "H={output.hist}) 2> {log}"
