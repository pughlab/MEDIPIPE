
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
You might come across above error even you spcified your input reads are pair-end. This might be due to your modified sample_template.tsv file is **not tab-delimited** or **with duplicated sample_id**. Fix that will allow the pipeline get throught. 


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

## "str' object has no attribute "name'
```bash
Traceback (most recent call last):
File "/cluster/home/arna/miniconda3/envs/MEDIPIPE/lib/python3.8/site-packages/snakemake/__init__py", line 699, in snakemake
success = workflow. execute(
File "/cluster/home/arna/miniconda/envs/MEDIPIPE/lib/python3.8/site-packages/snakemake/workflow.py", line 1056, in execute logger.run_info("\n" join(dag.statsO))
File "/cluster/home/arna/miniconda3/envs/MEDIPIPE/lib/python3.8/site-packages/snakemake/dag-py", line 2192, in stats yield tabulate(rows, headers="keys")
File "/cluster/home/arna/miniconda3/envs/MEDIPIPE/lib/python3.8/site-packages/tabulate/__init__.py", line 2048, in tabulate
list_of_lists, headers = _normalize_tabular_data
File "/cluster/home/arna/miniconda/envs/MEDIPIPE/lib/python3.8/site-packages/tabulate/__init__py", line 1471, in _normalize_tabular_data
rows = list(map(lambda r: r if _is_separating_line(r) else list(r), rows))
File "/cluster/home/arna/miniconda3/envs/MEDIPIPE/lib/python3.8/site-packages/tabulate/__init__-py", line 1471, in «lambda»
rows = list(map(lambda r: r if _is_separating line(r) else list(r), rows))
File "/cluster/home/arna/miniconda/envs/MEDIPIPE/lib/python3.8/site-packages/tabulate/__init__.py", line 107, in _is separating_line
(len (row) > 1 and row[0] = SEPARATING LINE)
File "/cluster/home/arna/miniconda3/envs/MEDIPIPE/lib/python3.8/site-packages/snakemake/rules.py", line 1127, in __eq-
return self.name = other.name and self.output = other.output
AttributeError: "str' object has no attribute "name'
```
This issure could be fix by "conda install tabulate=0.8.10", more details refer to snakemake issue [1892](https://github.com/snakemake/snakemake/issues/1892)

## Issue for installing mamba
You might encounter the issue when you intall the mamaba with your connda version. Downgrading your mamba to 1.#.# versions might fix it. The version used for building the pipeline is v0.15.3. 

## "ImportError: libffi.so.6: cannot open shared object file: No such file or directory" after OS upgrading 
This could be fixed by creating a softlink name that points to the exsiting file in the folder of /miniconda3/envs/MEDIPIPE/lib/,  more details refer to [here](https://stackoverflow.com/questions/61875869/ubuntu-20-04-upgrade-python-missing-libffi-so-6/63329830#63329830)
