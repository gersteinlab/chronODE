# Author: Mor Frank, with interface adaptations by Eve Wattenberg

# Imports
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset
from torch.utils.data import Dataset, DataLoader
import matplotlib.pyplot as plt
from torch.optim.lr_scheduler import LambdaLR
from sklearn.model_selection import train_test_split
import matplotlib
import scipy.stats as sc
import math
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputfile", type=str,
                help="Input file with paired RNA and OC values")
parser.add_argument("-p", "--predictfile", type=str,
                help="File to write predictions to")  
parser.add_argument("-S", "--subset", type=str, default="pos",
                help="Subset to model")                
parser.add_argument("-t", "--timepoints", type=str, default=105,
                help="Number of time points in the data")
parser.add_argument("-n", "--neuralnetfile", type=str, default=None,
                help="Filepath to save model, if used must end in .pth")                  
parser.add_argument("-m", "--modalities", type=str, default=1,
                help="Number of input data modalities")  
parser.add_argument("-s", "--split", type=str, default=0.2,
                help="Fraction of input to use as test set")                
parser.add_argument("-r", "--randomseed", type=str, default=1941,
                help="Random seed")        
args = parser.parse_args()

### Get data: Open chromatin signals and rna expression derivatives over time

# Read the databases of downregulated (decreasing) cCREs and upregulated (increasing) cCREs


# data_down = pd.read_csv('/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse/models/input/onetoone/downreg.filtered.merge.ESW.tsv',sep = '\t' )
# data_up = pd.read_csv('/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse/models/input/onetoone/upreg.filtered.merge.ESW.tsv',sep = '\t' )

# data = pd.concat([data_up,data_down], ignore_index=True)
data = pd.read_csv(args.inputfile, sep="\t")
data.drop_duplicates(inplace=True, keep='first', ignore_index = True)
# data.head()

## user needs to choose which model to run:
# positively correlated OC-rna Pairs (OC - open chromatin)
# negatively correalated 

print("which model do you want to run? for positively correlated OC-RNA type pos for negatively correlated type neg")
#user_input = input()
user_input = args.subset

print(f'user choice is {user_input}')

if user_input == 'pos':
    data_oc_rna = data[(data.correlation > 0)] 
if user_input == 'neg':
    data_oc_rna = data[(data.correlation < 0)]

#data_posR_posKrna
#data_posR_posKrna.head()
if user_input =='pos' or user_input =='neg':
    print(f'user chose {user_input}')
else:
    print('**** only pos or neg are allowed for model choices, please type again ****')

# data_oc_rna.head()

## Define: X_data (OC - open chromatin derivative features dc/dt) & 
## y_data (gene expression derivative target dg/dt)

def get_data(data):
    "input: dataframe with positive and negative correlations between oc and rna"
    "output: X_data, y_data for the model"

    #oc
    start_column_x = '10.5_oc'
    end_column_x = '21.0_oc'

    start_column_x = data.columns.get_loc(start_column_x)
    end_column_x = data.columns.get_loc(end_column_x)

    #rna
    start_column_y = '10.5_rna'      
    end_column_y = '21.0_rna'
    start_column_y = data.columns.get_loc(start_column_y)
    end_column_y = data.columns.get_loc(end_column_y)

    oc = data.iloc[:,start_column_x:end_column_x+1]
    rna = data.iloc[:,start_column_y:end_column_y+1]

    X_data = oc.values
    y_data = rna.values
    
    return (X_data, y_data)

X_data_OC_derivatives, y_data_rna_derivatives = get_data(data_oc_rna)

def tensorxy(X_data_OC_derivatives, y_data_rna_derivatives):
    "input: X (oc) and y(rna derivatives) arrays "
    "output: X (oc) and y(rna derivatives) tensors"  
    
    seq_len = 105 ## this is the number of time points


    X_OC_derivatives = X_data_OC_derivatives.reshape(len(X_data_OC_derivatives),1,seq_len)
    X_OC_derivatives = X_data_OC_derivatives.astype(float)


    seq_len_y = 105

    y_rna_derivatives = y_data_rna_derivatives.reshape(len(y_data_rna_derivatives),1,seq_len_y)

    # convert into PyTorch tensors
    
    #X = torch.tensor(X, dtype=torch.float32)
    X_OC_derivatives = torch.tensor(X_OC_derivatives, dtype=torch.float32)

    #y = torch.tensor(y, dtype=torch.float32)
    y_rna_derivatives = torch.tensor(y_rna_derivatives, dtype=torch.float32)
    
   
    return(X_OC_derivatives,y_rna_derivatives)

## Define X (Features) and y (target) as tensors
X1,y1 = tensorxy(X_data_OC_derivatives,y_data_rna_derivatives)
print("needles and pins")

### NN model

# TODO should this be parameterized?
batch_size_train = 4 ## can be changed to improve performance, larger--> runs faster. should be changed in train dataloader as well
batch_size_val = 1 
batch_size_test = 1 

#split the data into training set (80%) and temporary set (20%)


X_train1, X_test1, y_train1, y_test1, train_indexes1, test_indexes1 = train_test_split(X1, y1, np.arange(len(X1)),test_size=args.split, random_state=42)


n_points = args.timepoints

X_train1 = X_train1.reshape(-1, 1, n_points) 

X_test1 = X_test1.reshape(-1, 1, n_points) 

X_train1 = X_train1[:,0].reshape(X_train1.size()[0],1,n_points) # for OC feature only

train_dataset1 = TensorDataset(X_train1[:,:,:n_points], y_train1[:,:,:n_points])

train_loader1 = DataLoader(train_dataset1, batch_size=batch_size_train, shuffle=True)

X_test1 = X_test1[:,0].reshape(X_test1.size()[0],1,n_points) # for OC feature only

test_dataset1 = TensorDataset(X_test1[:,:,:n_points], y_test1[:,:,:n_points])

test_loader1 = DataLoader(test_dataset1, batch_size=batch_size_test, shuffle=False)


torch.manual_seed(42)
num_samples = y1.shape[0]


time_steps = args.timepoints
num_features = 1 # set number of features
batch_size = 4  # Define the batch size
loss_train_epoch = []

# Define network architecture to regress time series gene expression derivatives 
#from open chromatin signals
class PolynomialRegression(nn.Module):
    
    def __init__(self, degree, num_features):
        super(PolynomialRegression, self).__init__()
        self.degree = degree
        # set bias=True for posivitely correlated OC-rna derivative pairs
        # set bias=False for negatively correlated OC-rna derivative pairs 
        self.poly = nn.Linear(num_features, degree, bias=True)
        self.relu = nn.LeakyReLU(0.4,inplace=True)
        self.fc = nn.Linear(degree, 1)
 
        
    def forward(self, x):
        
        batch_size, time_steps, num_features = x.size()
        x_poly = self.poly(x)
        x_poly = x_poly.view(batch_size, time_steps, self.degree)
        x_poly = self.relu(x_poly)  # Apply Leaky ReLU activation
        out = self.fc(x_poly)
        
        return out

# Custom loss function to handle per-sample shape
class PerSampleMSELoss(nn.Module):
    def __init__(self):
        super(PerSampleMSELoss, self).__init__()
    

    def forward(self, input, target): 
        loss = torch.mean((input - target) ** 2, dim=1)  # Calculate the loss along the time_steps dimension
        
        return loss

# Initialize the model
degree = 30

warmup_steps = 1000  # Number of warm-up steps
initial_lr = 0.001  # Initial learning rate (warm-up learning rate)
desired_lr = 0.00001  # Desired learning rate

model = PolynomialRegression(degree, num_features)

# Define loss function and optimizer
criterion = PerSampleMSELoss()
optimizer = optim.Adam(model.parameters(), lr=initial_lr)

# Create a LambdaLR scheduler for warm-up
warmup_scheduler = LambdaLR(optimizer, lr_lambda=lambda step: step / warmup_steps if step < warmup_steps else 1)

# Training loop
num_epochs = 3000
for epoch in range(num_epochs):
 
    for i, (x_batch, y_batch) in enumerate(train_loader1):
        
        x_batch = x_batch.transpose(1, 2).requires_grad_()   ### mor #### NEED TO BE ADDED##########################
        y_batch = y_batch.transpose(1, 2).requires_grad_() 

        # Forward pass
        y_pred = model(x_batch)

        # Compute the loss
        loss = criterion(y_pred, y_batch)

        # Zero gradients, backward pass, and optimize
        optimizer.zero_grad()
        loss.mean().backward()
        optimizer.step()
        
        # Update the learning rate using the warm-up scheduler
        warmup_scheduler.step()

    if (epoch + 1) % 100 == 0:
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.mean().item():.4f}')
        loss_train_epoch.append(loss.mean().item())


## run this cell before loading a Pytorch model TODO removable?

## save model ###  ## to save a Pytorch model use:
#torch.save(model, 'model_name.pth')


class PolynomialRegression(nn.Module):
    def __init__(self, degree, num_features):
        super(PolynomialRegression, self).__init__()
        self.degree = degree
        # set bias=True for posivitely correlated OC-rna derivative pairs
        # set bias=False for negatively correlated OC-rna derivative pairs 
        self.poly = nn.Linear(num_features, degree, bias=False) 
        self.relu = nn.LeakyReLU(0.4,inplace=True)
        self.fc = nn.Linear(degree, 1)
     

    def forward(self, x):
        batch_size, time_steps, num_features = x.size()
        x_poly = self.poly(x)
        x_poly = x_poly.view(batch_size, time_steps, self.degree)
        x_poly = self.relu(x_poly)  # Apply Leaky ReLU activation
        
        
        out = self.fc(x_poly)
        
        return out
    
# Custom loss function to handle per-sample shape
class PerSampleMSELoss(nn.Module):
    def __init__(self):
        super(PerSampleMSELoss, self).__init__()
    
        

    def forward(self, input, target):
        loss = torch.mean((input - target) ** 2, dim=1)  # Calculate the loss along the time_steps dimension
        
        return loss
criterion = PerSampleMSELoss()
# TODO end of PyTorch block

# ### plot MSE (mean square error for each training epoch)

# plt.plot([i for i in range(0,3000,100)],loss_train_epoch)
# plt.xlabel('epoch')
# plt.ylabel('Avg. train MSE')

# Testing / inference
test_losses = []
all_pred_y = []
all_true_y = []

batch_size_test = 1
with torch.no_grad():
    
    for i, (x_batch, y_batch) in enumerate(test_loader1):
        
        x_batch = x_batch.transpose(1, 2)   
        y_batch = y_batch.transpose(1, 2)
        
        ## run this line if you trained the model
        #y_pred = model(x_batch)
        
        ## run this line to load the model for positively correlated OC-rna derivative pairs
        #model = torch.load('braindev/Model_Positive_OC_rna_corr.pth')
        
        ## run this line to load the model for negatively correlated OC-rna derivative pairs
        #model = torch.load('braindev/Model_Negative_OC_rna_corr.pth')
        
        model.eval()
        y_pred = model(x_batch)
        all_pred_y.append(y_pred)
        all_true_y.append(y_batch)

        # Compute the loss
       
        loss = criterion(y_pred, y_batch)
        print(f'Loss: {loss.mean().item():.4f}')
        test_losses.append(loss.mean().item())

print("All done! Saving!")

# save predictions
all_pred_y = np.squeeze(all_pred_y)
all_pred_y = pd.DataFrame(all_pred_y)
all_pred_y.to_csv(args.predictfile, sep="\t", index=False)
print("Saved predictions! Trying the model next")

## save model state
torch.save(model.state_dict(), 'model_name.pth')


# there were a huge number of plots here

# def find_indices(list_to_check, item_to_find):
#     "input: list_to_check: list if genes / cCREs ID"
#     "Output: item_to_find: gene / cCRE ID that its index needs to be found in list_to_check"
#     indices = []
#     for idx, value in enumerate(list_to_check):
#         if value == item_to_find:
#             indices.append(idx)
#     return indices

# # plot time series data of a specific cCRE along with its associated gene ID

# t = np.linspace(10.5,21,105)

# ## get list of all gene IDs which are in the text set
# l = [i[0] for i in data_oc_rna.iloc[test_indexes1[:],[1]].values.tolist()]

# ## find the index of a specific gene ID
# index_gene = find_indices(l, 'ENSMUSG00000015942')[0]
# print('gene index', find_indices(l, 'ENSMUSG00000015942'))
# geneid = f'{l[index_gene]}'
# print('gene id',l[index_gene])

# ## get list of all cCRE IDs which are in the text set 
# c = [i[0] for i in data_oc_rna.iloc[test_indexes1[:],[0]].values.tolist()]
# ccreid = f'{c[index_gene]}'
# print('cCRE id', c[index_gene])

# ## plot gene derivative true values versus predicted values
# plt.figure(figsize=(4,4))
# plt.plot(t, [i[0] for i in all_true_y[index_gene].tolist()[0]],'--ro')
# plt.plot(t, [i[0] for i in all_pred_y[index_gene].tolist()[0]], '--bo')

# ## print predicted gene derivatives
# print([i[0] for i in all_pred_y[index_gene][0].tolist()])


# print('first rna derivative value',data_oc_rna[(data_oc_rna['cCRE_id']==f'{c[index_gene]}') & (data_oc_rna['gene_id']==f'{l[index_gene]}')].loc[:,'10.5_rna'].tolist())
# data_oc_rna[(data_oc_rna['cCRE_id']==f'{c[index_gene]}') & (data_oc_rna['gene_id']==f'{l[index_gene]}')]

