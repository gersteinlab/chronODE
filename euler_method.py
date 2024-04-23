# Author: Mor Frank, with a new interface by Eve Wattenberg

# Imports
import pandas as pd
import numpy as np
import argparse
import sys

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

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--valuesfile", type=str,
                help="Input file with actual values")
parser.add_argument("-p", "--predictfile", type=str,
                help="Input file with predicted derivatives") 
parser.add_argument("-o", "--outfile", type=str,
                help="Output filepath")                
parser.add_argument("-t", "--timepoints", type=str, default=105,
                help="Number of time points in the data")  
parser.add_argument("-m", "--mode", type=str, default="all",
                help="Running mode")  
parser.add_argument("-g", "--geneid", type=str, default="",
                help="Gene ID to use with single gene mode")    
args = parser.parse_args()

print(args.mode)
if args.mode == "gene" and args.geneid == "":
    print("Please specify an ENSEMBL ID to use single-gene mode. Exiting.")
    sys.exit()

# Load "true" values
rna_actual_pipe = pd.read_csv(args.valuesfile,sep = '\t' )

if args.mode == "gene":
    real_values_geneid = args.geneid
    index_geneid = rna_actual_pipe[rna_actual_pipe['gene_id']== f'{real_values_geneid}'].index[0]
    initial_real_value_gene =  rna_actual_pipe.iloc[index_geneid,1]

# Load predicted derivatives
pred_data = pd.read_csv(args.predictfile, sep = '\t')
pred_data['gene_id'] = pred_data['combined_index'].str.split('-', expand=True)[0] #TODO parameterize this as optional???
pred_data.drop('combined_index', axis=1, inplace=True)

if args.mode == "gene":
    pred_index = pred_data[pred_data['gene_id']== f'{args.geneid}'].index[0]
    pred = pred_data.iloc[pred_index,0:-1].values.tolist()
else:
    # Select ONLY the genes that are present in both files
    intersect = rna_actual_pipe.merge(pred_data, on="gene_id")
    if len(intersect.index) == 0:
        print("No matching genes found. Check your file formats.")
        sys.exit()
    

# Construct line (gene mode)
if args.mode == "gene":
    constructed_line = run_euler(initial_real_value_gene, pred)

# Print output (gene mode)
if args.mode == "gene":
    print('true values', rna_actual_pipe[rna_actual_pipe['gene_id']== f'{real_values_geneid}'], '\n')
    print('reconstructed values', constructed_line[0:105], '\n')

# all mode 
if args.mode == "all":
    reconstructions = []
    for i in range(0, len(intersect.index)):
        initial_real_value_gene = intersect.iloc[i, 1]
        preds = intersect.iloc[i, (args.timepoints +1):]
        print("yes that slice worked...?")
        constructed_line = run_euler(initial_real_value_gene, preds)
        reconstructions.append(constructed_line)

    reconstructions = pd.DataFrame(reconstructions)
    print(reconstructions)
    print("Ciao!")

# Save results
if args.mode == "all":
    reconstructions.to_csv(args.outfile, sep="\t", index=False)
else:
    with open(args.outfile, 'w') as f:
        for val in constructed_line:
            f.write(f"{str(val)}\n")




## plot euler

# # TODO parameterize
# rna_actual_pipe = pd.read_csv('/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse/ODE/rna/values/forebrain.jan8.tsv',sep = '\t' )

# real_values_geneid = geneid 
# index_geneid = rna_actual_pipe[rna_actual_pipe['gene_id']== f'{real_values_geneid}'].index[0]
# real_values_gene =  rna_actual_pipe.iloc[index_geneid,1:]

#print('cCRE id', ccreid)
#print('gene id', real_values_geneid, '\n')


# plt.plot(t[0:105], constructed_line[0:105], linewidth=3, c='c')
# plt.plot(t[0:105],real_values_gene.tolist()[0:105], linewidth=3, c='r')

# plt.legend(['Reconstructed values', 'Real values'], fontsize=18, 
#            bbox_to_anchor=(0.64, 1.2), loc='center')
# plt.title(f'{ccreid}, {geneid}') 
# plt.xlabel('Time (PCDs)', fontsize=18)
# plt.ylabel('G', fontsize=18)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)

print("Goodbye!")