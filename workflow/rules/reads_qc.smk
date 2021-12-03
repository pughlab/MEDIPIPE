rule fastqc_before_trim:
    input:
        get_fastq
    output:
        "tmp/qc/{sample}_fastqc.html",
        "tmp/qc/{sample}_fastqc.zip"
    shell:
        "fastqc {input} --outdir tmp/qc/"

rule trim_galore_se:
    input:
        get_fastq
    output:
        temp("tmp/trimmed/{sample}_trimmed.fq"),                     ## suffix .fq.gz
        "tmp/trimmed/{sample}.fastq.gz_trimming_report.txt",
    params:
        path= config["workdir"] + "/tmp/trimmed/"                 ## path needs to be full path: failed to recognize space
    log:
        "tmp/trimmed/{sample}.log",
    shell:
        "trim_galore -q 20 --stringency 3 --length 20 --dont_gzip -o {params.path} {input}"             ## --dont_gzip is a walkaround of OneDirve - UNN

rule fastqc_after_trim:
    input:
        "tmp/trimmed/{sample}_trimmed.fq"
    output:
        "tmp/trimmed/{sample}_trimmed_fastqc.html",
        "tmp/trimmed/{sample}_trimmed_fastqc.zip"
    shell:
        "fastqc {input}"
