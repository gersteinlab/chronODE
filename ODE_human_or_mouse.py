### Author: ESW and Ke Xu
# DEPRECATED!!! USE ODE_fitting_posvals.py
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

# Function for fitting ODE to one time series
# if I understand this correctly, timepoints = np.linspace(start, end, num)
def subfit(vals, k_guess, b_guess, start, timepoints, positive):

    if positive:
        parameter_guesses = np.array([k_guess, b_guess])
    else:
        parameter_guesses = np.array([-k_guess, b_guess])
    
    def fitfunc(t, *params):

        def myode(x, t, parameters):
            k = parameters[0]
            b = parameters[1]
            dxdt = k * x * (1 - x/b)
            return dxdt

        # here, the vals should be the time series data for a certain cCRE/gene
        Xasol = odeint(myode, vals[0], t, args = (params,))
        return Xasol[:,0]

    try:
        fit_param, kcov = curve_fit(fitfunc, timepoints, vals, p0=parameter_guesses, maxfev = 5000)
    except Exception as e:
        print("curve_fit threw an exception for this input series.\nFilling in nans")
        y_ODE = np.full(len(timepoints), np.nan)
        y_ODE_cal = float("nan")
        fit_param = np.array([np.nan, np.nan])
        return y_ODE, y_ODE_cal, fit_param


    def system1order(y, t, k, b):
        dydt = k * y * (1 - y/b)
        return dydt

    fit = odeint(system1order, vals[0], timepoints, args=(fit_param[0], fit_param[1]))

    # find the analytical solution of the initial guess:
    k_ODE = fit_param[0]
    b_ODE = fit_param[1]
    c_ODE = math.log(abs(vals[0])+0.000001) - math.log(abs(1 - vals[0] / b_ODE)) - k_ODE * start
    C = np.e**c_ODE
    y_ODE_1 = (b_ODE * C * np.e**(k_ODE * timepoints))/(b_ODE + C * np.e**(k_ODE * timepoints))
    y_ODE_2 = - y_ODE_1

    # find mean squared error
    diff_ODE_1 = (y_ODE_1 - vals)**2
    diff_ODE_2 = (y_ODE_2 - vals)**2
    MSE_ODE_1 = diff_ODE_1.mean()
    MSE_ODE_2 = diff_ODE_2.mean()

    # choose between the two versions
    # the y_ODE_cal here will be used to define the openness in the later part!
    if MSE_ODE_1 <= MSE_ODE_2:
        y_ODE = y_ODE_1
        y_ODE_cal = (b_ODE * C * np.e**(k_ODE * (abs(1/k_ODE) + start)))/(b_ODE + C * np.e**(k_ODE * (abs(1/k_ODE) + start)))
    else:
        y_ODE = y_ODE_2
        y_ODE_cal = -(b_ODE * C * np.e**(k_ODE * (abs(1/k_ODE) + start)))/(b_ODE + C * np.e**(k_ODE * (abs(1/k_ODE) + start)))

    return y_ODE, y_ODE_cal, fit_param

def ODE_fit_sign(vals, k_guess, b_guess, start, timepoints):
        
    ########## dydt = -k * (1 - y / b) fit ##########
    y_ODE_neg, y_ODE_cal_neg, fit_param_neg = subfit(vals, k_guess, b_guess, start, timepoints, positive = False)

    ########## dydt = k * (1 - y / b) fit ###########
    y_ODE_pos, y_ODE_cal_pos, fit_param_pos  = subfit(vals, k_guess, b_guess, start, timepoints, positive = True)
        
        # Handle missing values. This is MESSY and might need tweaking
    if np.all(np.isnan(y_ODE_cal_neg)) and np.all(np.isnan(y_ODE_cal_pos)):
        print("both nan")
        # curve_fit failed twice so pass on the nans
        k_b_current = np.array([np.nan, np.nan]) 
        cal = np.nan
        T_cal = np.nan
        y_ODE = np.full(len(timepoints), np.nan)
        mse = np.nan
    elif np.all(np.isnan(y_ODE_cal_neg)):
        print("neg nan")
        # negative fit failed so just use the positive
        y_ODE = y_ODE_pos
        k_b_current = fit_param_pos
        cal = (y_ODE_cal_pos - y_ODE_pos.min()) / (y_ODE_pos.max() - y_ODE_pos.min())
        T_cal = abs(1/k_b_current[0])
        diff_guess_pos = (y_ODE_pos - vals)**2
        mse = diff_guess_pos.mean()
    elif np.all(np.isnan(y_ODE_cal_pos)):
        print("pos nan")
        # positive fit failed so just use negative
        y_ODE = y_ODE_neg
        k_b_current = fit_param_neg
        cal = (y_ODE_cal_neg - y_ODE_neg.min()) / (y_ODE_neg.max() - y_ODE_neg.min())
        T_cal = abs(1/k_b_current[0])
        diff_guess_neg = (y_ODE_neg - vals)**2
        mse = diff_guess_neg.mean()
    else:
        # define MSE between functional approach and the ODE approach

        diff_guess_neg = (y_ODE_neg - vals)**2
        diff_guess_pos = (y_ODE_pos - vals)**2

        MSE_guess_neg = diff_guess_neg.mean()
        MSE_guess_pos = diff_guess_pos.mean()

        if MSE_guess_neg < MSE_guess_pos:
            y_ODE = y_ODE_neg
            k_b_current = fit_param_neg
            # need to be very careful here when calculating the openness
            cal = (y_ODE_cal_neg - vals.min()) / (vals.max() - vals.min())
            T_cal = abs(1/k_b_current[0])
            mse = MSE_guess_neg
        else:
            y_ODE = y_ODE_pos
            k_b_current = fit_param_pos
            cal = (y_ODE_cal_pos - vals.min()) / (vals.max() - vals.min())
            T_cal = abs(1/k_b_current[0])
            mse = MSE_guess_pos
                
    return(k_b_current, y_ODE, mse, cal, T_cal)

# function with move, this is the main fitting function
def ODE_fit(df, num_timepoints, t_orig):
    # Defining some constants
    start = t_orig[0]
    end = t_orig[-1]
    timepoints = np.linspace(start, end, num_timepoints)
    k_guess = 0.9
    b_guess = 1.5
    
    # Build a dataframe of interpolated values for each gene, also do the min_max normalization here
    df_normalized = df.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=1)
    interp_vals = []
    # the following is based on the min-max normalization version
    for i in range(0,len(df_normalized)):
        interp_func = interp1d(t_orig, df_normalized.iloc[i,:].tolist())
        interp_down_sig = interp_func(timepoints)
        interp_vals.append(interp_down_sig)
        
    ### fitting the parameters ###
    k_b = [] # fitted parameters
    calculation = [] # to be discussed
    T = [] # to be discussed
    MSE = []
    TYPE = [] # up, down, org
    MOVE = []
    fitted_vals = [] # the fitted values for the fitted function
    fitted_derivs = [] # don't sure if we still need this
    
    for i in range(0,len(interp_vals)):
        current_vals = interp_vals[i] # This corresponds to the red dots in Ke's graphs
        move = np.max(abs(current_vals))
        vals_upward = current_vals + move
        vals_downward = current_vals - move
        
        outcome_org = ODE_fit_sign(current_vals, k_guess, b_guess, start, timepoints)
        outcome_up = ODE_fit_sign(vals_upward, k_guess, b_guess, start, timepoints)
        outcome_down = ODE_fit_sign(vals_downward, k_guess, b_guess, start, timepoints)

        # List of outcomes with respective type labels and adjustment values
        outcomes = [
            (outcome_up, 'upward', -move),
            (outcome_down, 'downward', move),
            (outcome_org, 'original', 0)
        ]
        
        # Selecting the outcome with the smallest MSE
        best_outcome, Type, adjustment = min(outcomes, key=lambda x: x[0][2])
        k_b_best, fit, mse, cal, T_cal = best_outcome
        
        # Adjust the fit
        fit += adjustment
        # calculate the derivates for each fit
        derivs = k_b_best[0] * fit * (1 - fit / k_b_best[1]) # Ke's version
        
        # Add values to list
        k_b.append(k_b_best)
        calculation.append(cal)
        T.append(T_cal)
        MSE.append(mse)
        TYPE.append(Type)
        MOVE.append(move)
        fitted_vals.append(fit) # This corresponds to the yellow lines in Ke's graphs
        fitted_derivs.append(derivs)
    
    # Assemble parameters DF
    k_b = pd.DataFrame(k_b, columns = ['k', 'b'])
    prop = pd.DataFrame(calculation, columns=['proportion'])
    Time = pd.DataFrame(T, columns=['T'])
    MSE = pd.DataFrame(MSE, columns=['MSE'])
    TYPE = pd.DataFrame(TYPE, columns=['TYPE'])
    MOVE = pd.DataFrame(MOVE, columns=['MOVE'])
    index = df.index
    k_b = k_b.set_axis(index, axis=0)
    prop = prop.set_axis(index, axis=0)
    Time = Time.set_axis(index, axis=0)
    MSE = MSE.set_axis(index, axis=0)
    TYPE = TYPE.set_axis(index, axis=0)
    MOVE = MOVE.set_axis(index, axis=0)
    parameters = pd.concat([k_b, prop, Time, MSE,TYPE,MOVE], axis = 1)

    # Assemble modeled values DF
    fitted_vals = pd.DataFrame(fitted_vals)
    fitted_vals.columns = timepoints
    fitted_vals.index = df.index

    # Assemble modeled derivatives DF
    fitted_derivs = pd.DataFrame(fitted_derivs)
    fitted_derivs.columns = timepoints
    fitted_derivs.index = df.index
    
    return(parameters, fitted_vals, fitted_derivs, timepoints)
        

if __name__ == "__main__":
    print("starting")
    os.chdir("/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse/ODE/bin/")
    print("parsing")
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
                    help="assay suffix to put on column names")
    parser.add_argument("-v", "--valuesfile", type=str,
                    help="Values output file")
    parser.add_argument("-d", "--derivsfile", type=str,
                    help="Derivatives output file")
    parser.add_argument("-p", "--paramsfile", type=str,
                    help="Parameters output file")
    args = parser.parse_args()

    print(f"Loading {args.assay} data from {args.group} {args.region}, using file\n{args.inputfile}")
    print(f"Modeling values and derivatives for {args.timepoints} interpolated timepoints.")
    print(f"Values save to {args.valuesfile}\nDerivatives save to {args.derivsfile}\nODE parameters save to {args.paramsfile}")

    # Load data
    raw_data = pd.read_csv(args.inputfile, sep = '\t', index_col=0)

    # Select timecourse
    if args.timecourse == "mouse":
        original_timecourse = [10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 21]
    else if args.timecourse == "human":
        original_timecourse = [56, 72, 80, 85, 96, 104, 105, 109, 112, 117, 122] # I hope this is correct
    else: 
        print("Timecourse " + args.timecourse + " not supported. Please choose mouse or human")
        sys.exit()

    # run fitting
    parameters, values, derivatives, interp_timepoints = ODE_fit(raw_data, args.timepoints)

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

    values.to_csv(args.valuesfile, sep='\t')
    derivatives.to_csv(args.derivsfile, sep = '\t')
    parameters.to_csv(args.paramsfile, sep='\t')

    print("Successfully saved to")
    print(args.derivsfile)