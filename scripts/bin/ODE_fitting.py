#!/usr/bin/env python

# Author: ESW and Ke
# BB made some minor edits
# This version of the code requires y > 0 and b > 0

import pandas as pd
from scipy import signal
from scipy.interpolate import interp1d
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import odeint
import math
import sys
import os
import argparse

# This is the core fitting function
def subfit(myvals, k_guess, b_guess, start, timepoints, positive):

    # Decide whether to test k>0 or k<0
    if positive:
        parameter_guesses = np.array([k_guess, b_guess])
    else:
        parameter_guesses = np.array([-k_guess, b_guess])

    # Define ODE function
    def fitfunc(t, *params):
        def myode(y, t, parameters):
            k = parameters[0]
            b = parameters[1]
            dydt = k * y * (1 - y/b)
            return dydt

        Xasol = odeint(myode, myvals[0], t, args=(params,), mxstep = 5000)
        return Xasol[:,0]

    # Estimate k and b
    try:
        fit_param, kcov = curve_fit(fitfunc, timepoints, myvals, p0=parameter_guesses, maxfev=5000)
    except Exception as e:
        print("curve_fit threw an exception for this input series.\nFilling in nans")
        return np.full(len(timepoints), np.nan), np.nan, [np.nan, np.nan]

    # Use estimated k and b to compute C
    k_ODE = fit_param[0]
    b_ODE = fit_param[1]

    try:
        c_ODE = math.log(abs(myvals[0]) + 0.000001) - math.log(abs(1 - myvals[0] / b_ODE)) - k_ODE * start
        C = np.exp(c_ODE)

        # Find the analytical solution (length = timepoints)
        y_ODE = (b_ODE * C * np.exp(k_ODE * timepoints)) / (b_ODE + C * np.exp(k_ODE * timepoints))

        # Compute mean squared error
        diff_ODE = (y_ODE - myvals) ** 2
        MSE_ODE = np.nanmean(diff_ODE)

        # Return the solution if b > 0
        if b_ODE > 0:
            return y_ODE, MSE_ODE, fit_param 
        else:
            return np.nan, np.nan, [np.nan, np.nan]
    except:
        return np.full(len(timepoints), np.nan), np.nan, [np.nan, np.nan]


        
# This function is a wrapper for subfit
# It tries both positive and negative k and
# selects the fitting with lowest MSE
def ODE_fit_sign(input_vals, k_guess, b_guess, start, timepoints):

    # fit using k < 0
    y_ODE_neg, MSE_ODE_neg, fit_param_neg = subfit(input_vals, k_guess, b_guess, start, timepoints, positive=False)

    # fit using k > 0 
    y_ODE_pos, MSE_ODE_pos, fit_param_pos = subfit(input_vals, k_guess, b_guess, start, timepoints, positive=True)

    # pick the most successful fit (i.e. the one that is not NA and has the lowest MSE)
    if np.isnan(MSE_ODE_neg) and not np.isnan(MSE_ODE_pos):
        return fit_param_pos, y_ODE_pos, MSE_ODE_pos 
    elif np.isnan(MSE_ODE_pos) and not np.isnan(MSE_ODE_neg):
        return fit_param_neg, y_ODE_neg, MSE_ODE_neg 
    elif not np.isnan(MSE_ODE_pos) and not np.isnan(MSE_ODE_neg):
        if MSE_ODE_neg <= MSE_ODE_pos:
            return fit_param_neg, y_ODE_neg, MSE_ODE_neg 
        else:
            return fit_param_pos, y_ODE_pos, MSE_ODE_pos
    else:
        return np.nan, np.nan, np.nan 


# This function is a wrapper for ODE_fit_sign
# It runs ODE_fit_sign for every element in the input dataframe
def ODE_fit(df, t_orig):

    # Define start & end of the timecourse
    start = t_orig[0]
    end = t_orig[-1]

    # Define initial guesses for k and b
    k_guess = 0.9
    b_guess = 1.5

    # Perform fitting
    vals = [row.to_numpy(dtype=float) for _, row in df.iterrows()]
    results = []
    for val in vals:  
        
        # Fit the ODE
        outcome = ODE_fit_sign(val, k_guess, b_guess, start, t_orig) 

        # Save fitting results; if unsuccessful fitting, return NA
        if not np.isnan(outcome[2]):
            results.append((outcome))
        else:
            results.append((np.nan, np.nan, np.nan, np.nan))

    # Unpack results into a DataFrame
    parameters = pd.DataFrame({
        'k': [result[0][0] if not np.isnan(result[2]) else np.nan for result in results],
        'b_starred': [result[0][1] if not np.isnan(result[2]) else np.nan for result in results],
        'MSE': [result[2] for result in results]
     }, index=df.index)

    # Create fitted_vals DataFrame
    fitted_vals = pd.DataFrame({
        tp: [r[1][i] if isinstance(r[1], np.ndarray) and not np.isnan(r[1][i]) else np.nan for r in results]
        for i, tp in enumerate(t_orig)
    }, index=df.index)

    # Calculate derivatives for each fitting
    fitted_derivs = pd.DataFrame({
        tp: [
            result[0][0] * val[i] * (1 - val[i] / result[0][1]) 
            if isinstance(val, np.ndarray) and not np.isnan(val[i]) 
            else np.nan 
            for result, val in zip(results, fitted_vals.values)
        ] 
        for i, tp in enumerate(t_orig)
    }, index=df.index)

    return parameters, fitted_vals, fitted_derivs

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--inputfile", type=str,
                    help="Input file with values at 8 timepoints")

    parser.add_argument("-t", "--timepoints", type=str,
                    help="Set of initial timepoints to use")

    parser.add_argument("-v", "--valuesfile", type=str,
                    help="Values output file")

    parser.add_argument("-d", "--derivsfile", type=str,
                    help="Derivatives output file")

    parser.add_argument("-p", "--paramsfile", type=str,
                    help="Parameters output file")

    args = parser.parse_args()

    # Load input data
    input_data = pd.read_csv(args.inputfile, sep = '\t', index_col=0)

    # Read in timecourse
    original_timecourse = pd.read_csv(args.timepoints, header=None)
    original_timecourse = np.array(original_timecourse.values.flatten().tolist())

    # Run fitting
    parameters, values, derivatives = ODE_fit(input_data, original_timecourse)

    # Change column names for clarity
    def colnames_format(t):
        return(str(round(t, 3)))
    colnames = list(map(colnames_format, original_timecourse))
    values.columns = colnames
    derivatives.columns = colnames

    # Replace nans with "NA" for readable output
    values = values.fillna("NA")
    derivatives = derivatives.fillna("NA")
    parameters = parameters.fillna("NA")

    values.to_csv(args.valuesfile, sep='\t')
    derivatives.to_csv(args.derivsfile, sep = '\t')
    parameters.to_csv(args.paramsfile, sep='\t')
