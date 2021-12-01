rule bwa_mem_pe:
    input:
        reads=["reads/{sample}.1.fastq", "reads/{sample}.2.fastq"]
    output:
        "mapped/{sample}.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index="genome",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="none",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "0.71.0/bio/bwa/mem"
        
        
rule bwa_mem_se:
    input:
        reads=["reads/{sample}.1.fastq"]
    output:
        "mapped/{sample}.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index="genome",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="none",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "0.71.0/bio/bwa/mem"
        
        
rule samtools_sort:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    params:
        extra = "-m 4G",
        tmp_dir = "/tmp/"
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "0.71.0/bio/samtools/sort"
        
rule samtools_index:
    input:
        "mapped/{sample}.sorted.bam"
    output:
        "mapped/{sample}.sorted.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.71.0/bio/samtools/index"
        
