
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
