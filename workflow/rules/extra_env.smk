################################################################################
## Check and install extra conda envs automatically,
## Need internet during the inital installation and updating

### all extra env
rule install_all_extra_env:
    input:
        'extra_env/bedtools_installed',
        'extra_env/ConsensusCruncher_installed',
        'extra_env/multiQC_installed',
        'extra_env/R_installed',
        'extra_env/umi_tools_installed'
    output:
        'extra_env/all_extra_env_installed'
    shell:
        'touch {output}'

## individual extra env
rule install_extra_env_4_bedtools:
    output:
        'extra_env/bedtools_installed'
    conda:
        'extra_env/bedtools.yaml'
    shell:
        'touch {output}'

rule install_extra_env_4_ConsensusCruncher:
    output:
        'extra_env/ConsensusCruncher_installed'
    conda:
        'extra_env/ConsensusCruncher.yaml'
    shell:
        'touch {output}'

rule install_extra_env_4_multiQC:
    output:
        'extra_env/multiQC_installed'
    conda:
        'extra_env/multiQC.yaml'
    shell:
        'touch {output}'

rule install_extra_env_4_R:
    output:
        'extra_env/R_installed'
    conda:
        'extra_env/R.yaml'
    shell:
        'touch {output}'

rule install_extra_env_4_umi_tools:
    output:
        'extra_env/umi_tools_installed'
    conda:
        'extra_env/umi_tools.yaml'
    shell:
        'touch {output}'
