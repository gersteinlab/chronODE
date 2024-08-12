#!/usr/bin/env python

# Adapted by ESW from BB's jupyter notebook

# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from sklearn.metrics import mean_squared_error
import argparse

def permute(f_vals, o_vals, permutations=10000) :

    # compute original MSE
    original_mse = mean_squared_error(f_vals, o_vals)

    # compute permuted MSE
    shuffled_mses = []

    for i in range(permutations):
        np.random.seed(i)
        shuffled_vals = np.random.permutation(o_vals)
        shuffled_mse = mean_squared_error(f_vals, shuffled_vals)
        shuffled_mses.append(shuffled_mse)

    p_value = np.mean(np.array(shuffled_mses) <= original_mse)
    return(p_value)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-r", "--rescvals", type=str,
                    help="Input file with rescaled fitted values")

    parser.add_argument("-o", "--origvals", type=str,
                    help="Input file with original values")

    parser.add_argument("-p", "--paramsfile", type=str,
                    help="Input file with fitted parameters")

    args = parser.parse_args()


# Load params 
m_params = pd.read_csv(args.paramsfile, sep="\t")
m_params_table = np.array(m_params)

# Load original values
m_origin_vals = pd.read_csv(args.origvals, sep="\t")
m_origin_vals_table = np.array(m_origin_vals.iloc[:,1:m_origin_vals.shape[1]])

# Load rescaled fitted values
m_fitted_vals = pd.read_csv(args.rescvals, sep="\t")
m_fitted_vals_table = np.array(m_fitted_vals.iloc[:,1:m_fitted_vals.shape[1]])

# Loop through elements and compute permutation p-value for each element
new_params = []

for i in range(len(m_params_table)):
    params = m_params_table[i].tolist()

    o_vals = m_origin_vals_table[i].tolist()
    f_vals = m_fitted_vals_table[i].tolist()

    pval = permute(f_vals, o_vals)

    params.append(pval)
    new_params.append(params)


column_names = m_params.columns.tolist()
column_names.append('p_value')
new_params = pd.DataFrame(new_params, columns = column_names)

# Save results
new_params.to_csv("mse.param.tsv", sep="\t", index=False)




