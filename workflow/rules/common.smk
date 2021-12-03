
def get_fastq(wildcards):
    return config["samples"][wildcards.sample]
