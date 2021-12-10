    layout = SAMPLES.loc[wildcards.sample]["layout"],

    if layout == "se":
        R1 = SAMPLES.loc[wildcards.sample]["R1"],
        return R1
    else:
        R1 = SAMPLES.loc[wildcards.sample]["R1"],
        R2 = SAMPLES.loc[wildcards.sample]["R2"],
        return R1 + R2


rule fastqc_before_trim:
    input:
        get_fastq
    output:
        "tmp/qc/{sample}_fastqc.html",
        "tmp/qc/{sample}_fastqc.zip"
    shell:
        "fastqc {input} --outdir tmp/qc/"


rule fastqc_after_trim:
    input:
        "tmp/trimmed/{sample}_trimmed.fq"
    output:
        "tmp/trimmed/{sample}_trimmed_fastqc.html",
        "tmp/trimmed/{sample}_trimmed_fastqc.zip"
    shell:
        "fastqc {input}"





"""
    if config["paired-end"]:
        tmp1 = expand("trimmed_fq/{sample}_R1_val_1.fq", sample = SAMPLES["sample"]),
        tmp2 = expand("trimmed_fq/{sample}_R2_val_2.fq", sample = SAMPLES["sample"]),
        return tmp1 + tmp2
    else:
        return expand("trimmed_fq/{sample}_trimmed.fq", sample = SAMPLES["sample"]),
"""
