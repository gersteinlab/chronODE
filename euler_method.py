# Author: Mor Frank, with a new interface by Eve Wattenberg

# Imports
import pandas as pd
print("i am running")

def euler(t_prev,y_prev, deriv_prev,t_next):
    "apply the Euler method"
    "input: t_prev: previous time point (t_n) ,y_prev: previous gene / cCRE value"
    "deriv_prev: derivative at time t_n, t_next: the next time point (t_n+1)"
    
    y_next = y_prev + deriv_prev*(t_next - t_prev)
    
    return y_next

def run_euler(initial_point, all_y_pred):
    "run the Euler method"
    "input: initial point: real value of gene expression"
    "all_y_pred: gene expression at all time points"
     
    t = np.linspace(10.5,21,105)
    n = 0
    constructed_line = []
    y=initial_point  ## initial point is the real point
    y0=[initial_point]

    for i in range(0,105):
      
        y_pred = all_y_pred[n]
        t_prev = t[n]
        y_prev = y
        deriv_prev = y_pred
        if n+1==105:
            break
        t_next = t[n+1]
        y = euler(t_prev,y_prev, deriv_prev, t_next)
        constructed_line.append(y)
        n = n+1
    constructed_line = y0+constructed_line 
    return constructed_line

# get parameters to run the Euler method

print("like a good little script")
geneid = "ENSMUSG00000015942"
#cCRE id EM10D1960585

# TODO parameterize
rna_actual_pipe = pd.read_csv('/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse/ODE/rna/values/forebrain.jan8.tsv',sep = '\t' )
real_values_geneid = geneid
index_geneid = rna_actual_pipe[rna_actual_pipe['gene_id']== f'{real_values_geneid}'].index[0]
initial_real_value_gene =  rna_actual_pipe.iloc[index_geneid,1]

### run euler for the NN

# TODO paameterize
rna_actual_pipe = pd.read_csv('/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse/ODE/rna/values/forebrain.jan8.tsv',sep = '\t' )
real_values_geneid = geneid
index_geneid = rna_actual_pipe[rna_actual_pipe['gene_id']== f'{real_values_geneid}'].index[0]
initial_real_value_gene =  rna_actual_pipe.iloc[index_geneid,1]

all_y_pred = all_pred_y[index_gene].tolist()[0] 
all_y_pred = [i[0]for i in all_y_pred]
constructed_line = run_euler(initial_real_value_gene, all_y_pred)

# ## run euler for the RF

# rna_actual_pipe = pd.read_csv('/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse/ODE/rna/values/forebrain.jan8.tsv',sep = '\t' )
# real_values_geneid = geneid
# index_geneid = rna_actual_pipe[rna_actual_pipe['gene_id']== f'{real_values_geneid}'].index[0]
# initial_real_value_gene =  rna_actual_pipe.iloc[index_geneid,1]

# RF_pred_data = pd.read_csv('rf_predictions_samples.tsv', sep = '\t')
# #RF_pred = RF_pred_data.iloc[4,1:-1].tolist()  ### put the row number
# RF_pred = RF_pred_data[RF_pred_data['gene_id'] == geneid].iloc[0,1:-1].values.tolist()
# constructed_line = run_euler(initial_real_value_gene, RF_pred)


## plot euler

# TODO parameterize
rna_actual_pipe = pd.read_csv('/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse/ODE/rna/values/forebrain.jan8.tsv',sep = '\t' )

real_values_geneid = geneid 
index_geneid = rna_actual_pipe[rna_actual_pipe['gene_id']== f'{real_values_geneid}'].index[0]
real_values_gene =  rna_actual_pipe.iloc[index_geneid,1:]

print('cCRE id', ccreid)
print('gene id', real_values_geneid, '\n')
print('true values', real_values_gene.tolist(), '\n')
print('reconstructed values', constructed_line[0:105], '\n')

# plt.plot(t[0:105], constructed_line[0:105], linewidth=3, c='c')
# plt.plot(t[0:105],real_values_gene.tolist()[0:105], linewidth=3, c='r')

# plt.legend(['Reconstructed values', 'Real values'], fontsize=18, 
#            bbox_to_anchor=(0.64, 1.2), loc='center')
# plt.title(f'{ccreid}, {geneid}') 
# plt.xlabel('Time (PCDs)', fontsize=18)
# plt.ylabel('G', fontsize=18)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)


### NN performance: Pearson's r an MAE (mean absolute error)

## plot true-predicted correlation histogram

corr_list = []
true_list = []
pred_list = []

for i in range(0, len(all_true_y)):
    true_list.append([i[0] for i in all_true_y[i][0].tolist()])
    pred_list.append([i[0] for i in all_pred_y[i][0].tolist()])

for i in range(0, len(pred_list)):
    r,p = sc.pearsonr(pred_list[i], true_list[i])
    corr_list.append(r)

corr_list = [0 if math.isnan(x) else x for x in corr_list]
# plt.figure(figsize=(5,5))
# plt.hist(corr_list, histtype = 'stepfilled', color='royalblue')
# plt.xticks(fontsize=18)
# plt.xlabel('true - predicted dg/dt \n correlation', fontsize=18)
# plt.yticks(fontsize=18)
# plt.ylabel('counts', fontsize=18)
# plt.ylim(0,350)
# #plt.xlim(-0.1,1)
# plt.xticks([-1,-0.5, 0, 0.5, 1])
# #plt.yticks([0, 100])
print(np.nanpercentile(corr_list, 75))
print(np.nanpercentile(corr_list, 50))
print(np.nanpercentile(corr_list, 25)) 

print(np.mean(corr_list))
print(np.std(corr_list))
print('positive corr', len([i for i in corr_list if i>0])/len(corr_list))