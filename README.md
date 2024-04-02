# *chronODE*: A framework to integrate time-series multi-omics data based on ordinary differential equations combined with machine learning
Beatrice Borsari, Mor Frank, Eve S. Wattenberg, Susanna X. Liu, Xuezhu Yu, Mark Gerstein  
[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.12.13.571513v1)  

***

## Requirements

    pandas=2.2.1
    

## Example running
Run pipeline on open chromatin data from mouse forebrain decreasing cCREs:  

    python ODE_fitting_posvals.py \
      --inputfile "example_dir/open_chromatin_example_input.tsv" \
      --timepoints 105 \
      --timecourse "mouse" \
      --group "decreasing" \
      --region "forebrain" \
      --assay "open_chromatin" \
      --valuesfile "example_dir/oc_vals.tsv" \
      --derivsfile "example_dir/oc_derivs.tsv" \
      --paramsfile "example_dir/oc_params.tsv"
Run pipeline on DNase-seq data from human brain:  

    python ODE_fitting_posvals.py \
      --inputfile "example_dir/human_rnaseq.tsv" \
      --timepoints 105 \
      --timecourse "human" \
      --assay "DNase" \
      --valuesfile "example_dir/dnase_vals.tsv" \
      --derivsfile "example_dir/dnase_derivs.tsv" \
      --paramsfile "example_dir/dnase_params.tsv"
