
## automatically remove adapters for single-end reads
rule trim_galore_se:
    input:
        "reads/{sample}.fastq.gz"
    output:
        "trimmed/{sample}_trimmed.fq.gz",
         "trimmed/{sample}.fastq.gz_trimming_report.txt"
    params:
        extra="--illumina -q 20"
    log:
        "logs/trim_galore/{sample}.log"
    wrapper:
        "0.71.0/bio/trim_galore/se"

## automatically remove adapters for paired-end reads
rule trim_galore_pe:
    input:
        ["reads/{sample}.1.fastq.gz", "reads/{sample}.2.fastq.gz"]
    output:
        "trimmed/{sample}.1_val_1.fq.gz",
         "trimmed/{sample}.1.fastq.gz_trimming_report.txt",
         "trimmed/{sample}.2_val_2.fq.gz",
         "trimmed/{sample}.2.fastq.gz_trimming_report.txt"
    params:
        extra="--illumina -q 20"
    log:
        "logs/trim_galore/{sample}.log"
    wrapper:
        "0.71.0/bio/trim_galore/pe"

## Fastqc report for reads after removing adapters
rule fastqc:
    input:
        unpack(get_fastq),
    output:
        html="results/qc/fastqc/{sample}-{unit}.html",
        zip="results/qc/fastqc/{sample}-{unit}.zip",
    log:
        "logs/fastqc/{sample}-{unit}.log",
    wrapper:
        "0.74.0/bio/fastqc"

