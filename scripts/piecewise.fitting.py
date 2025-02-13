## import packages
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
from sklearn.metrics import mean_squared_error
import torch
import multiprocessing as mp
import os
from sklearn.metrics import auc


####   Computing peak-like fitting    ########


def min_max(vals, c, d):
	"min-max normalization"  
    
    c = c
    d = d
    min_vals = min(vals)
    max_vals = max(vals)
    vals = ((d-c)*(vals - min_vals)) / (max_vals - min_vals) + c
    
    return vals


def find_local_maxmin(vals):
    "input: all signal values"
    "output: 2 signals. left side and right side with the corresponded time points"
    
    
    timepoints = [10.5, 11.5, 12.5, 13.5, 14.5,15.5, 16.5, 21]

    
    max_vals = np.argmax(vals)   
    min_vals = np.argmin(vals)     

    avg_max_min = (vals[max_vals]+vals[min_vals])/2
    
   
    norm_array=vals  

        left_side = []
        right_side =[]
        max_vals = []
        min_vals = []
        flag = 1
        
    else:
    
        if vals[max_vals] - np.mean(vals[max_vals]-vals[min_vals]) > np.mean(vals) - vals[min_vals]:

            left_side = vals[0:(max_vals+1)]
            right_side = vals[(max_vals):]

            flag = 1

        else:                     ## opens upward
            
            if min_vals < len(vals):
                left_side = vals[0:(min_vals+1)]
            else:
                left_side = vals[0:(min_vals)]

            right_side = vals[(min_vals):]

            flag = 0
    
    return timepoints, left_side, right_side, max_vals, min_vals, flag


def fit_left_right(timepoints,left_side, right_side, max_vals, min_vals):
	"piecewise fitting"

    if flag == 1:

        vals_left = left_side
        vals_left = vals_left
        try:
            timepoints_left = timepoints[0:(max_vals+1)]
        except:
            timepoints_left = []

        vals_right = right_side
        vals_right = vals_right
        try:
            timepoints_right = timepoints[max_vals:]
        except:
            timepoints_right = []

    else:

        vals_left = left_side
        vals_left = vals_left
        
        try:
            timepoints_left = timepoints[0:min_vals+1]
        except:
            timepoints_left = []
            
        vals_right = right_side
        vals_right = vals_right
        try:
            timepoints_right = timepoints[min_vals:] 
        except:
            timepoints_right = []
        
    return timepoints_left, timepoints_right, vals_left, vals_right  



# define fit function: this will find k & b
def fitfunc(t, *params):
        def myode(y, t, parameters):
            k = parameters[0]
            b = parameters[1]
            dydt = k * y * (1 - y/b) 
            return dydt

        Xasol = odeint(myode, vals[0], t, args=(params,), mxstep=5000, printmessg=False)
        return Xasol[:,0]


## define function to find the fitted curve 
## using k & b parameters found from previous code
# def system1order(y, t, k, b):
#     dydt = k*y*(1- y/b) + t
#     return dydt

def system1order(y, t, k, b):
    dydt = k * y * (1 - y/b) 
    
    return dydt


def perfrom_fitting(time_points, values, system1order, fit_param):
    "input: left or right timepoints and values"
    "output: fitted values for side of the curve"
    
    # this step returns a vector of length = 8 with the values of the fitted curve
    fit_side = odeint(system1order, values[0], time_points, args=(fit_param[0], fit_param[1]), printmessg=False)

    return time_points, fit_side


filename = 'File_Name'   ### Peak file name
m = pd.read_csv(f'{filename}.tsv', sep='\t')

m.head()

cCREid = m.iloc[:,0].values.tolist()
m.shape[0]


# ### run the cell below to find fitted lines , k and b in the 0-1 and 1-2 ranges

number_signals = m.shape[0]
with mp.Pool(processes=30) as pool:
    for c,d in [[0.00001, 1], [1,2]]:   

        timepoints = [10.5, 11.5, 12.5, 13.5, 14.5,15.5, 16.5, 21]

        b2guess_b = 1.5
        succeed_left = 0

        save_fit = []
        save_fit_left = []
        save_fit_right = []
        save_time_left = []
        save_time_right = []
        save_true_sig = []
        save_param_right = []
        save_param_left = []


        ids = []

        my_array = np.array(m.iloc[:,1:9])

        ## min max normalization between 0 to 1
        for i in range(0, number_signals):                   

            norm_array = min_max(my_array[i], c, d)

            ## find minima / maxima
            timepoints, left_side, right_side, max_vals, min_vals, flag = find_local_maxmin(norm_array)

            ## divide the signal to left and right: time point and signal values
            timepoints_left, timepoints_right, vals_left, vals_right  = fit_left_right(timepoints,left_side, right_side, max_vals, min_vals)

            if flag == 1:
                b2guess_k = 0.9

            else:
                b2guess_k = -0.9

            ## make a fitting for each side: left / right
            #if flag == 1:
            "fit left side k and b"
            parameterguesses = np.array([b2guess_k , b2guess_b])
            vals = vals_left

            try:
                fit_param_left, kcov = curve_fit(fitfunc, timepoints_left, vals, p0=parameterguesses , maxfev=5000)
                succeed_left = 1
                save_param_left.append(fit_param_left)

            except:

                succeed_left = 0
                save_fit_left.append([])
                save_time_left.append([])
                save_fit_right.append([])
                save_time_right.append([])

                save_param_right.append([])
                save_param_left.append([])

                #save_true_sig.append(norm_array)
                save_true_sig.append(min_max(my_array[i], c, d))
                ids.append(cCREid[i])
                continue

            timepoints_left, fit_left = perfrom_fitting(timepoints_left, vals_left, system1order, fit_param_left)
            save_fit_left.append(fit_left.reshape(1,len(fit_left))[0].tolist())
            save_time_left.append(timepoints_left)


            if succeed_left == 1:

                "fit right side k and b"
                parameterguesses = np.array([-b2guess_k, b2guess_b])
                vals = vals_right

                vals[0] = fit_left.reshape(1,len(fit_left))[0].tolist()[-1] ## remove discontinuity


                try:
                    fit_param_right, kcov = curve_fit(fitfunc, timepoints_right, vals, p0=parameterguesses , maxfev=5000)
                    save_param_right.append(fit_param_right)
                except:

                    succeed_left = 0
                    save_fit_right.append([])
                    save_time_right.append([])

                    save_param_right.append([])

                    save_true_sig.append(norm_array)
                    ids.append(cCREid[i])
                    continue

                timepoints_right, fit_right = perfrom_fitting(timepoints_right, vals_right, system1order, fit_param_right)
                save_fit_right.append(fit_right.reshape(1,len(fit_right))[0].tolist())
                save_time_right.append(timepoints_right)


            else:
                succeed_left = 0
                save_fit_right.append([])
                save_time_right.append([])

                save_param_right.append([])

            save_true_sig.append(min_max(my_array[i], c, d))
            ids.append(cCREid[i])

            print(f'******** signal {i}, flag {flag}, fitted parameters left {fit_param_left}, fitted parameters right {fit_param_right}********')
        if d==1:    
            save_time_left01 = save_time_left
            save_fit_left01 = save_fit_left
            save_time_right01 = save_time_right
            save_fit_right01 = save_fit_right
            timepoints01 = timepoints
            save_true_sig01 = save_true_sig
            save_param_right01 = save_param_right
            save_param_left01 = save_param_left
    pool.close()
    pool.join()


### Compute MSE
def MSE_calc(fitted_left, fitted_right, true_sig, param_left, param_right):
	"MSE computation"
    
    cum_MSE = []
    cum_fitted_lines = []
    
    
    for i in range(0, len(true_sig)): 
        try:
            fitted_line = fitted_left[i][:-1]+fitted_right[i]
            if param_left[i][1]>0 and param_right[i][1]>0 and np.sign(param_left[i][0])!=np.sign(param_right[i][0]):
                cum_fitted_lines.append(fitted_line)
                cum_MSE.append(mean_squared_error(true_sig[i], fitted_line))
            else:
                cum_MSE.append([])
                cum_fitted_lines.append([])
        
        except:
            cum_MSE.append([])
            cum_fitted_lines.append([])
            continue
                           
    return cum_MSE, cum_fitted_lines


MSE01, fitted_line01 = MSE_calc(save_fit_left01, save_fit_right01, save_true_sig01, save_param_left01, save_param_right01)
MSE12, fitted_line12 = MSE_calc(save_fit_left, save_fit_right, save_true_sig, save_param_left, save_param_right)


# ## save best fitting 0-1 or 1-2 range

## save the best fit 

save_best_fit = []
true_lines = []
best_koc_left = []
best_koc_right = []

for i in range(0, len(MSE01)):
    
    if MSE12[i]==[]:
        

        true_lines.append(save_true_sig01[i])

        if save_param_left01[i]!=[] and save_param_right01[i]!=[]:
            if save_param_left01[i][1]>0 and save_param_right01[i][1]>0:

                save_best_fit.append(fitted_line01[i])
                best_koc_left.append(save_param_left01[i])
                best_koc_right.append(save_param_right01[i])

            else:
                save_best_fit.append([])
                best_koc_left.append([])
                best_koc_right.append([])
        else:
            save_best_fit.append([])
            best_koc_left.append([])
            best_koc_right.append([])

        
    elif MSE01[i]==[]:


        true_lines.append(save_true_sig[i])
        if save_param_left[i]!=[] and  save_param_right[i]!=[]:
            if save_param_left[i][1]>0 and save_param_right[i][1]>0:
            
                save_best_fit.append(fitted_line12[i])
                best_koc_left.append(save_param_left[i])
                best_koc_right.append(save_param_right[i])
            else:
                save_best_fit.append([])
                best_koc_left.append([])
                best_koc_right.append([])
        else:
            save_best_fit.append([])
            best_koc_left.append([])
            best_koc_right.append([])


    elif MSE01[i] < MSE12[i]:

            true_lines.append(save_true_sig01[i])
            if save_param_right01[i]!=[] and save_param_left01[i]!=[]:
                if save_param_right01[i][1]>0 and save_param_left01[i][1]>0:
                    save_best_fit.append(fitted_line01[i])
                    best_koc_left.append(save_param_left01[i])
                    best_koc_right.append(save_param_right01[i])
                else:
                    save_best_fit.append([])
                    best_koc_left.append([])
                    best_koc_right.append([])
            else:
                save_best_fit.append([])
                best_koc_left.append([])
                best_koc_right.append([])

    else:
            true_lines.append(save_true_sig[i])
            if save_param_left[i]!=[] and save_param_right[i]!=[]:
                if save_param_left[i][1]>0 and save_param_right[i][1]>0:
                    save_best_fit.append(fitted_line12[i])
                    best_koc_left.append(save_param_left[i])
                    best_koc_right.append(save_param_right[i])
                else:
                    save_best_fit.append([])
                    best_koc_left.append([])
                    best_koc_right.append([])
            else:
                save_best_fit.append([])
                best_koc_left.append([])
                best_koc_right.append([])


            try:
                if (best_koc_left[i][1] < save_true_sig[i][-1]) or (best_koc_right[i][1] < save_best_fit[i][0]) and flag==1:
                    save_best_fit[-1] = []
                    best_koc_left[-1] = []
                    best_koc_right[-1] = []
            except:
                if (best_koc_left[i][1] < save_true_sig[i][0]) or (best_koc_right[i][1] < save_best_fit[i][-1]) and flag==0:
                    save_best_fit[-1] = []
                    best_koc_left[-1] = []
                    best_koc_right[-1] = []

# ### save best fitting in dataframe

timepoints = [10.5, 11.5, 12.5, 13.5, 14.5,15.5, 16.5, 21]

def save_fittings_in_dataframe(ids, true_lines, save_best_fit):
	"save normalized signals"

    d = {'gene_id': ids, 
         '10.5_true': [i[0] for i in true_lines], 
         '11.5_true': [i[1] for i in true_lines], 
         '12.5_true': [i[2] for i in true_lines], 
         '13.5_true': [i[3] for i in true_lines], 
         '14.5_true': [i[4] for i in true_lines], 
         '15.5_true': [i[5] for i in true_lines], 
         '16.5_true': [i[6] for i in true_lines], 
         '21_true': [i[7] for i in true_lines],
         
         '10.5_fitted': [i[0] if i!=[] else [] for i in save_best_fit], 
         '11.5_fitted':[i[1] if i!=[] else [] for i in save_best_fit], 
         '12.5_fitted': [i[2] if i!=[] else [] for i in save_best_fit], 
         '13.5_fitted':[i[3] if i!=[] else [] for i in save_best_fit], 
         '14.5_fitted': [i[4] if i!=[] else [] for i in save_best_fit], 
         '15.5_fitted':[i[5] if i!=[] else [] for i in save_best_fit], 
         '16.5_fitted':[i[6] if i!=[] else [] for i in save_best_fit], 
         '21_fitted':[i[7] if i!=[] else [] for i in save_best_fit],
         
         
         'k_rna_left': [i[0] if i!=[] else [] for i in best_koc_left],
         'k_rna_right':[i[0] if i!=[] else [] for i in best_koc_right],
         'b_rna_left': [i[1] if i!=[] else [] for i in best_koc_left],
         'b_rna_right':[i[1] if i!=[] else [] for i in best_koc_right]}
    return pd.DataFrame(data=d)


# ### save dataframe to csv file

## save to dataframe

my_peaks = save_fittings_in_dataframe(ids, true_lines, save_best_fit)

# # Save SAVE normalized table (0-1 / 1-2)
my_peaks.to_csv('normalaized_table.csv')
my_peaks.to_csv(f'{filename}_extreme_p_norm.csv')


# ### read csv between 0-1 / 1-2 before doing rescaling
m_fitted = pd.read_csv(f'{filename}_extreme_p_norm.csv', index_col=0)


# ## save the table in rescaled format (reverse from 0-1 / 1-2 to the original range)

### table with original values

fitted_rescaled = []


for i in range(0, len(m_fitted)):  #len(m_fitted)
    
    vals_norm = m_fitted.iloc[i, 9:17]
    try:
        
        x = list(map(float, vals_norm))
        minval = m.iloc[i, 1:9].min()
        if np.max(x) > 1.11:
            maxval = (m.iloc[i,1:9].max() -m.iloc[i,1:9].min()) * (np.max(x)-1) + m.iloc[i,1:9].min()
            minval = (m.iloc[i,1:9].max() -m.iloc[i,1:9].min()) * (np.min(x)-1) + m.iloc[i,1:9].min()

        else:
            maxval =  (m.iloc[i,1:9].max() -m.iloc[i,1:9].min()) * (np.max(x)) + m.iloc[i,1:9].min()
            minval = (m.iloc[i,1:9].max() -m.iloc[i,1:9].min()) * (np.min(x)) + m.iloc[i,1:9].min()
        fitted_rescaled.append(min_max(np.array(x), minval, maxval))
    
    except:
        fitted_rescaled.append([])


#replace columns in dataframe to the rescaled (no normalization)
rescaled_array = fitted_rescaled

col0 = pd.Series([i[0] if i!=[] else [] for i in rescaled_array])
col1 = pd.Series([i[1] if i!=[] else [] for i in rescaled_array])
col2 = pd.Series([i[2] if i!=[] else [] for i in rescaled_array])
col3 = pd.Series([i[3] if i!=[] else [] for i in rescaled_array])
col4 = pd.Series([i[4] if i!=[] else [] for i in rescaled_array])
col5 = pd.Series([i[5] if i!=[] else [] for i in rescaled_array])
col6 = pd.Series([i[6] if i!=[] else [] for i in rescaled_array])
col7 = pd.Series([i[7] if i!=[] else [] for i in rescaled_array])


col_names = m.columns[1:].tolist()
m_fitted['10.5_true'] = m[col_names[0]]
m_fitted['11.5_true'] = m[col_names[1]]
m_fitted['12.5_true'] = m[col_names[2]]
m_fitted['13.5_true'] = m[col_names[3]]
m_fitted['14.5_true'] = m[col_names[4]]
m_fitted['15.5_true'] = m[col_names[5]]
m_fitted['16.5_true'] = m[col_names[6]]
m_fitted['21_true'] = m[col_names[7]]
                       
m_fitted['10.5_fitted'] = col0
m_fitted['11.5_fitted'] = col1
m_fitted['12.5_fitted'] = col2
m_fitted['13.5_fitted'] = col3
m_fitted['14.5_fitted'] = col4
m_fitted['15.5_fitted'] = col5
m_fitted['16.5_fitted'] = col6
m_fitted['21_fitted'] = col7                          


# ### SAVE table rescaled (not normalized)

m_fitted.to_csv(f'{filename}_extreme_p_rescale.csv')

# read again the normalized table
m_fitted = pd.read_csv(f'{filename}_extreme_p_norm.csv', index_col = 0)

# ### save final table
rescaled_table.to_csv(f'{filename}_extreme_p_rescale.csv')
m_fitted.to_csv(f'{filename}_extreme_p_norm.csv')