# *chronODE*: A framework to integrate time-series multi-omics data based on ordinary differential equations combined with machine learning
Beatrice Borsari, Mor Frank, Eve S. Wattenberg, Susanna X. Liu, Xuezhu Yu, Mark Gerstein  
[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.12.13.571513v1)  

An ODE-based  pipeline for interpolating values and derivatives from time-series data, and machine learning models that predicted gene expression from linked chromatin features.
A figure should go here!

***

## *chronODE*

### Requirements

```Libraries should go here?```
### Input data and files format

```
Usage: ODE_fitting_posvals.py [-h] [-i INPUTFILE] [-t TIMEPOINTS] [-T TIMECOURSE] [-g GROUP] [-r REGION] [-a ASSAY] [-v VALUESFILE] [-d DERIVSFILE]
                              [-p PARAMSFILE] [-f FOLDER]
Required arguments:
  -i INPUTFILE, --inputfile INPUTFILE
                        Input file with values at a handful of timepoints
  -v VALUESFILE, --valuesfile VALUESFILE
                        Values output file
  -d DERIVSFILE, --derivsfile DERIVSFILE
                        Derivatives output file
  -p PARAMSFILE, --paramsfile PARAMSFILE
                        Parameters output file
  -f FOLDER, --folder FOLDER
                        Base filepath for file names
optional arguments:
  -h, --help            show this help message and exit

  -t TIMEPOINTS, --timepoints TIMEPOINTS
                        Number of timepoints to interpolate; defaults to 105
  -T TIMECOURSE, --timecourse TIMECOURSE
                        Set of initial timepoints to use, either mouse (default) or human
  -g GROUP, --group GROUP
                        expression pattern; defauts to 'unknown'
  -r REGION, --region REGION
                        Mouse brain region; defaults to 'unknown'
  -a ASSAY, --assay ASSAY
                        Assay suffix to put on column names; defaults to 'x'

```

### Examples
Run pipeline on RNAseq data from mouse forebrain:  

    python ODE_fitting_posvals.py \
      --inputfile "example_dir/rna_example_input.tsv" \
      --timepoints 105 \
      --timecourse "mouse" \
      --region "forebrain" \
      --assay "RNAseq" \
      --valuesfile "example_dir/rna_vals.tsv" \
      --derivsfile "example_dir/rna_derivs.tsv" \
      --paramsfile "example_dir/rna_params.tsv"
Run pipeline on DNase-seq data from human brain decreasing cCREs:  

    python ODE_fitting_posvals.py \
      --inputfile "example_dir/human_dnase.tsv" \
      --timepoints 105 \
      --timecourse "human" \
      --group "decreasing" \
      --assay "DNase" \
      --valuesfile "example_dir/dnase_vals.tsv" \
      --derivsfile "example_dir/dnase_derivs.tsv" \
      --paramsfile "example_dir/dnase_params.tsv"
