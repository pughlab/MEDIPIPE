
## Conda environment
It might pop out the compile error as below when you activate the MEDIPIPE, it is might be due to base conda env. However, you can still successfully run the pipeline for now. Will try to fix this soon. 

```bash
ERROR: This cross-compiler package contains no program /cluster/home/yzeng/miniconda3/env
s/MEDIPIPE/bin/x86_64-conda_cos6-linux-gnu-gfortran
INFO: activate-gfortran_linux-64.sh made the following environmental changes:
+HOST=x86_64-conda_cos6-linux-gnu
-HOST=x86_64-conda-linux-gnu
```

## Missing input files for rule merge_and_rename_fq_se
You might come across above error even you spcified your input reads are pair-end. This might be due to your modified sample_template.tsv file is not tab-delimited or with duplicated sample_id. Fix that will allow the pipeline get throught. 


## Could not create conda environment from umi_tools.yaml
```bash
LibMambaUnsatisfiableError: Encountered problems while solving:
  - nothing provides libgcc-ng >=12 needed by samtools-1.17-h00cdaf9_0
 
Could not solve for environment specs
The following package could not be installed
└─ samtools 1.17**  is not installable because it requires
   └─ libgcc-ng >=12 , which does not exist (perhaps a missing channel).
```
This issue could be fixed by adding channels "conda-forge" and "defaults" into umi_tools.yaml, which has also been updated in the latest version.
