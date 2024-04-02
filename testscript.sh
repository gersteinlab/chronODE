#!/bin/bash
#SBATCH --job-name=example
#SBATCH --output="test_april2.out"
#SBATCH --time=1:00:00
#SBATCH --partition=pi_gerstein
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL

# Author: ESW

# load in python environments (you may not need this step)
module load miniconda
conda activate jupyter_env

# WARNING: This path is based on our local Yale filepath!
example_dir="/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/paper/github/chronODE/example_data"

# Run chronODE pipeline on open chromatin
python ODE_fitting_posvals.py \
  --inputfile "$example_dir/open_chromatin_example_input.tsv" \
  --timepoints 105 \
  --timecourse "mouse" \
  --group "decreasing" \
  --region "forebrain" \
  --assay "open_chromatin" \
  --valuesfile "$example_dir/oc_vals.tsv" \
  --derivsfile "$example_dir/oc_derivs.tsv" \
  --paramsfile "$example_dir/oc_params.tsv"

# Run chronODE on RNAseq
python ODE_fitting_posvals.py \
  --inputfile "$example_dir/rnaseq_example_input.tsv" \
  --timepoints 105 \
  --timecourse "mouse" \
  --region "forebrain" \
  --assay "rnaseq" \
  --valuesfile "$example_dir/rnaseq_vals.tsv" \
  --derivsfile "$example_dir/rnaseq_derivs.tsv" \
  --paramsfile "$example_dir/rnaseq_params.tsv"