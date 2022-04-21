###############################################
## get a copy of raw fq and rename to smaple ID
## to make sure consistant wildcard.sample
################################################
## single end
rule rename_fq_se:
    input:
        get_raw_fastq
    output:
        temp("renamed_fq/{sample}.fastq.gz"),
    shell:
        "cp {input} {output} "

## paired-end
rule rename_fq_pe:
    input:
        get_raw_fastq
    output:
        temp("renamed_fq/{sample}_R1.fastq.gz"),
        temp("renamed_fq/{sample}_R2.fastq.gz"),
    shell:
        "cp {input[0]} {output[0]} && "
        "cp {input[1]} {output[1]} "


###########################################################
### extract UMI barcode and add it to FASTQ headers, p.r.n.
### pair-end unzipped FASTQ only, so far
###########################################################
rule extract_barcode:
    input:
        get_renamed_fastq
    output:
        temp("barcoded_fq/{sample}_R1.fastq"),
        temp("barcoded_fq/{sample}_R2.fastq"),
        "barcoded_fq/{sample}_R1.fastq.gz",
        "barcoded_fq/{sample}_R2.fastq.gz"
    conda:
        "extra_env/ConsensusCruncher.yaml"
    params:
        src = pipe_dir + "/workflow/dependencies/ConsensusCruncher/ConsensusCruncher/extract_barcodes.py",
        blist = umi_list,                 ## barcode list
        outfile = "barcoded_fq/{sample}"
    shell:
        ## unzip gz files
        "gunzip {input[0]} -c > {params.outfile}_R1.fastq && "
        "gunzip {input[1]} -c > {params.outfile}_R2.fastq && "

        ## extract barcodes
        #"python  {params.src} --bpattern NNT "
        "python  {params.src} --blist {params.blist} "
        "--read1  {input[0]}  --read2   {input[1]} "
        "--outfile  {params.outfile} && "

        ## gzip
        "gzip  {params.outfile}_barcode_R*.fastq  && "

        ## change to consistant names for following steps
        "mv {params.outfile}_barcode_R1.fastq.gz  {params.outfile}_R1.fastq.gz && "
        "mv {params.outfile}_barcode_R2.fastq.gz  {params.outfile}_R2.fastq.gz"



###################################
### automatically trimming adapters
### -q 20  Quality gz_trimming
### mimimal length : 20 nt
###################################
#for single-end reads
rule trim_galore_se:
    input:
        get_fastq_4trim
    output:
        temp("trimmed_fq/{sample}_trimmed.fq.gz"),
        "trimmed_fq/{sample}.fastq.gz_trimming_report.txt",
    params:
        ## path needs to be full path: failed to recognize space
        path = wd + "/trimmed_fq"
    threads: 12
    log:
        "logs/{sample}_trim_galore_se.log"
    shell:
        "(trim_galore -q 20 --stringency 3 --length 20 "
        "--cores {threads} -o {params.path} {input}) 2> {log}"

#for paired-end reads
rule trim_galore_pe:
    input:
        get_fastq_4trim
    output:
        temp("trimmed_fq/{sample}_R1_val_1.fq.gz"),
        temp("trimmed_fq/{sample}_R2_val_2.fq.gz"),
        "trimmed_fq/{sample}_R1.fastq.gz_trimming_report.txt",
        "trimmed_fq/{sample}_R2.fastq.gz_trimming_report.txt",
    params:
        ## path needs to be full path: failed to recognize space
        path = wd + "/trimmed_fq"
    threads: 12
    log:
        "logs/{sample}_trim_galore_pe.log"
    shell:
        "(trim_galore -q 20 --stringency 3 --length 20 "
        "--cores {threads} --paired -o {params.path} {input}) 2> {log}"



##########################################
### FASTQC for raw and trimmed fastq reads
##########################################
#for single-end reads
rule fastqc_se:
    input:
        get_renamed_fastq,
        get_trimmed_fastq
    output:
        "fastqc_se/{sample}_fastqc.zip",
        "fastqc_se/{sample}_trimmed_fastqc.zip"
    run:
        for fq in input:
            shell("fastqc {} -t 8 --outdir fastqc_se/".format(fq))

#for paired-end reads
rule fastqc_pe:
    input:
        get_renamed_fastq,
        get_trimmed_fastq
    output:
        "fastqc_pe/{sample}_R1_fastqc.zip",
        "fastqc_pe/{sample}_R2_fastqc.zip",
        "fastqc_pe/{sample}_R1_val_1_fastqc.zip",
        "fastqc_pe/{sample}_R2_val_2_fastqc.zip"
    run:
        for fq in input:
            shell("fastqc {} -t 8 --outdir fastqc_pe/".format(fq))
