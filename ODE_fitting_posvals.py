# Authors: Ke Xu and Eve Wattenberg
# This version requires y > 0 and b > 0

# Imports
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
print("imported")

# Core ODE fitting function
def subfit(vals, k_guess, b_guess, start, timepoints, positive):
    # Choose initial parameter set
    if positive:
        parameter_guesses = np.array([k_guess, b_guess])
    else:
        parameter_guesses = np.array([-k_guess, b_guess])

    # Define the ODE to fit
    def fitfunc(t, *params):
        def myode(x, t, parameters):
            k = parameters[0]
            b = parameters[1]
            dxdt = k * x * (1 - x/b)
            return dxdt

        Xasol = odeint(myode, vals[0], t, args=(params,))
        return Xasol[:,0]

    try:
        fit_param, kcov = curve_fit(fitfunc, timepoints, vals, p0=parameter_guesses, maxfev=5000)
    except Exception as e:
        print("curve_fit threw an exception for this input series.\nFilling in nans")
        return np.full(len(timepoints), np.nan), np.nan, np.nan, [np.nan, np.nan] 

    # def system1order(y, t, k, b):
    #     dydt = k * y * (1 - y/b)
    #     return dydt

    # find the analytical solution of the initial guess:
    k_ODE = fit_param[0]
    b_ODE = fit_param[1]
    c_ODE = math.log(abs(vals[0]) + 0.000001) - math.log(abs(1 - vals[0] / b_ODE)) - k_ODE * start
    C = np.exp(c_ODE)
    y_ODE_p = (b_ODE * C * np.exp(k_ODE * timepoints)) / (b_ODE + C * np.exp(k_ODE * timepoints))
    y_ODE_n = -y_ODE_p

    # find mean squared error
    diff_ODE_p = (y_ODE_p - vals) ** 2
    diff_ODE_n = (y_ODE_n - vals) ** 2
    MSE_ODE_p = np.nanmean(diff_ODE_p)
    MSE_ODE_n = np.nanmean(diff_ODE_n)

    # choose between the two versions [FLAGGING CHANGE] old code decides this differently
    if b_ODE > 0:
        if MSE_ODE_p <= MSE_ODE_n:
            return y_ODE_p, MSE_ODE_p, 1, fit_param
        else:
            return y_ODE_n, MSE_ODE_n, -1, fit_param
    else:
        return np.nan, np.nan, np.nan, [np.nan, np.nan]

# Wrapper function that tries fits with negative
# and positive k and selects the best fit
def ODE_fit_sign(vals, k_guess, b_guess, start, timepoints):
    # negative initial k
    y_ODE_neg, MSE_ODE_neg, sign_func_neg, fit_param_neg = subfit(vals, k_guess, b_guess, start, timepoints, positive=False)
    # positive initial k
    y_ODE_pos, MSE_ODE_pos, sign_func_pos, fit_param_pos = subfit(vals, k_guess, b_guess, start, timepoints, positive=True)

    # pick the only successful fit
    if np.isnan(MSE_ODE_neg) and not np.isnan(MSE_ODE_pos):
        return fit_param_pos, y_ODE_pos, MSE_ODE_pos, sign_func_pos
    elif np.isnan(MSE_ODE_pos) and not np.isnan(MSE_ODE_neg):
        return fit_param_neg, y_ODE_neg, MSE_ODE_neg, sign_func_neg
    # if both fits succeeded, select the one with lowest MSE 
    elif not np.isnan(MSE_ODE_pos) and not np.isnan(MSE_ODE_neg):
        if MSE_ODE_neg <= MSE_ODE_pos:
            return fit_param_neg, y_ODE_neg, MSE_ODE_neg, sign_func_neg
        else:
            return fit_param_pos, y_ODE_pos, MSE_ODE_pos, sign_func_pos
    # if both fits failed, pass on nans
    else:
        return np.nan, np.nan, np.nan, np.nan

# Main fitting function
def ODE_fit(df, num_timepoints, t_orig): 
    # Define the interpolated timecourse
    start = t_orig[0]
    end = t_orig[-1]
    timepoints = np.linspace(start, end, num_timepoints)

    # Initial guesses for k and b
    k_guess = 0.9
    b_guess = 1.5

    # Min-max normalization of input
    df_normalized = df.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=1)

    # Interpolate input values to match timecourse
    interp_vals = [interp1d(t_orig, row, bounds_error=False, fill_value='extrapolate')(timepoints) for row in df_normalized.values] 

    # Perform fitting!
    results = []
    for vals in interp_vals: 
        # Create a version of the input values that is guaranteed to be all positive
        move = np.max(abs(vals))
        vals_upward = vals + move

        # Fit ODE to original and shifted input values
        outcome_org = ODE_fit_sign(vals, k_guess, b_guess, start, timepoints)
        outcome_up = ODE_fit_sign(vals_upward, k_guess, b_guess, start, timepoints)

        # Select the best result of the two versions:

        # pick the only successful fit
        if np.isnan(outcome_up[2]) and not np.isnan(outcome_org[2]):
            results.append((*outcome_org, 'original', 0))
        elif np.isnan(outcome_org[2]) and not np.isnan(outcome_up[2]):
            results.append((*outcome_up, 'upward', move))
        # if both fits succeeded, select the one with lowest MSE
        elif not np.isnan(outcome_up[2]) and not np.isnan(outcome_org[2]):
            if outcome_up[2] <= outcome_org[2]: 
                results.append((*outcome_up, 'upward', move))
            else:
                results.append((*outcome_org, 'original', 0))
        # if both fits failed, pass on nans
        else:
            results.append((np.nan, np.nan, np.nan, np.nan, 'NA', np.nan))

    # Unpack results into a DataFrame
    parameters = pd.DataFrame({
        'k': [result[0][0] if not np.isnan(result[2]) else np.nan for result in results],
        'b': [result[0][1] if not np.isnan(result[2]) else np.nan for result in results],
        'MSE': [result[2] for result in results],
        'sign_func': [result[3] for result in results],
        'TYPE': [result[4] for result in results],
        'MOVE': [result[5] for result in results]
    }, index=df.index)

    # Create fitted values DataFrame
    fitted_vals = pd.DataFrame({
        tp: [r[1][i] if isinstance(r[1], np.ndarray) and not np.isnan(r[1][i]) else np.nan for r in results]
        for i, tp in enumerate(timepoints)
    }, index=df.index)

    # Calculate derivatives for each fitting
    fitted_derivs = pd.DataFrame({
        tp: [
            result[3] * result[0][0] * val[i] * (1 - val[i] / result[0][1]) 
            if isinstance(val, np.ndarray) and not np.isnan(val[i]) 
            else np.nan 
            for result, val in zip(results, fitted_vals.values)
        ] 
        for i, tp in enumerate(timepoints)
    }, index=df.index)

    return parameters, fitted_vals, fitted_derivs, timepoints

# Main body
if __name__ == "__main__":
    print("starting")
    print("parsing")

    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputfile", type=str,
                    help="Input file with values at 8 timepoints")
    parser.add_argument("-t", "--timepoints", type=int,
                    help="Number of timepoints to interpolate")
    parser.add_argument("-T", "--timecourse", type=str, default="mouse",
                    help="Set of initial timepoints to use")
    parser.add_argument("-g", "--group", type=str, default="unknown",
                    help="expression pattern")
    parser.add_argument("-r", "--region", type=str, default="unknown",
                    help="Mouse brain region")
    parser.add_argument("-a", "--assay", type=str, default="x",
                    help="Assay suffix to put on column names")
    parser.add_argument("-v", "--valuesfile", type=str,
                    help="Values output file")
    parser.add_argument("-d", "--derivsfile", type=str,
                    help="Derivatives output file")
    parser.add_argument("-p", "--paramsfile", type=str,
                    help="Parameters output file")
    parser.add_argument("-f", "--folder", type=str, 
                    default="/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse/ODE/",
                    help="Base folder")
    args = parser.parse_args()

    # Sanity-check messages to ensure you have used the correct CLI options
    print(f"Loading {args.assay} data from {args.group} {args.region}, using file\n{args.inputfile}")
    print(f"Modeling values and derivatives for {args.timepoints} interpolated timepoints.")
    print(f"Values save to {args.valuesfile}\nDerivatives save to {args.derivsfile}\nODE parameters save to {args.paramsfile}")

    # Load raw data
    os.chdir(args.folder)
    raw_data = pd.read_csv(args.inputfile, sep = '\t', index_col=0)

    # Select timecourse
    if args.timecourse == "mouse":
        original_timecourse = [10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 21]
    elif args.timecourse == "human":
        original_timecourse = [56, 72, 80, 85, 96, 104, 105, 109, 112, 117, 122, 142] 
    else: 
        print("Timecourse " + args.timecourse + " not supported. Please choose mouse or human")
        sys.exit()

    # run fitting
    parameters, values, derivatives, interp_timepoints = ODE_fit(raw_data, args.timepoints, original_timecourse)

    parameters['group'] = args.group
    parameters['region'] = args.region

    # Change column names for clarity
    def colnames_format(t):
        t = str(round(t, 3))
        return(t + '_' + args.assay)
    colnames = list(map(colnames_format, interp_timepoints))
    values.columns = colnames
    derivatives.columns = colnames

    # Replace nans with "NA" for readable output
    values = values.fillna("NA")
    derivatives = derivatives.fillna("NA")
    parameters = parameters.fillna("NA")

    # Save values, derivatives, and parameters
    values.to_csv(args.valuesfile, sep='\t')
    derivatives.to_csv(args.derivsfile, sep = '\t')
    parameters.to_csv(args.paramsfile, sep='\t')

    print("Successfully saved to")
    print(args.derivsfile)