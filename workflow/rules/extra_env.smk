################################################################################
## Check and install extra conda envs automatically,
## Need internet during the inital installation and updating
## wildcard {sample} : installed

### all extra env
rule install_all_extra_env:
    input:
        'extra_env/bedtools_{sample}',
        'extra_env/ConsensusCruncher_{sample}',
        'extra_env/multiQC_{sample}',
        'extra_env/R_{sample}',
        'extra_env/r_aggr_qc_report_{sample}',
        'extra_env/umi_tools_{sample}'
    output:
        'extra_env/all_extra_env_{sample}'    ## using wildcard works for --cluster as well
    shell:
        'touch {output}'

## individual extra env
rule install_extra_env_4_bedtools:
    output:
        'extra_env/bedtools_{sample}'
    conda:
        'extra_env/bedtools.yaml'
    shell:
        'touch {output}'

rule install_extra_env_4_ConsensusCruncher:
    output:
        'extra_env/ConsensusCruncher_{sample}'
    conda:
        'extra_env/ConsensusCruncher.yaml'
    shell:
        'touch {output}'

rule install_extra_env_4_multiQC:
    output:
        'extra_env/multiQC_{sample}'
    conda:
        'extra_env/multiQC.yaml'
    shell:
        'touch {output}'

rule install_extra_env_4_R:
    output:
        'extra_env/R_{sample}'
    conda:
        'extra_env/R.yaml'
    shell:
        'touch {output}'

rule install_extra_env_4_R_aggr_qc_report:
    output:
        'extra_env/r_aggr_qc_report_{sample}'
    conda:
        'extra_env/R_aggr_qc_report.yaml'
    shell:
        'touch {output}'


rule install_extra_env_4_umi_tools:
    output:
        'extra_env/umi_tools_{sample}'
    conda:
        'extra_env/umi_tools.yaml'
    shell:
        'touch {output}'
