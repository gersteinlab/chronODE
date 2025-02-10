#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset
from torch.utils.data import Dataset, DataLoader
import matplotlib.pyplot as plt
from torch.optim.lr_scheduler import LambdaLR
from sklearn.model_selection import train_test_split
import scipy.stats as sc
import math
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import iqr
from torch import nn
import torch.optim as optim


# In[3]:


torch.cuda.manual_seed_all(42)
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
device


# ### Get data: Open chromatin signals and rna expression over time

# In[4]:


data = pd.read_csv('****.tsv', sep = '\t')
data.head()


# In[5]:


## number of cCREs per time points
num_cCREs_per_time_p = 60     ###60


# In[7]:


def save_multi_image(filename):
    pp = PdfPages(filename)
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()


# In[10]:


## min max normalization 

def min_max_norm(my_array):
    
    normalized_data = []

    for i in my_array:
        
            try:
                min_val = np.min(i, axis =1).reshape(num_cCREs_per_time_p,1)
                max_val = np.max(i, axis =1).reshape(num_cCREs_per_time_p,1)
            except:

                min_val = np.min(i)
                max_val = np.max(i)

            
            norm_vals = (i- min_val) / (max_val - (min_val))
            arr_cleaned = np.nan_to_num(norm_vals, nan=0.0)
            normalized_data.append(arr_cleaned)
            
        
    return normalized_data


# In[9]:


### Choose regulatory mechanism


#data_oc_rna = data[(data['group'] == 'group3_repressor') | (data['group'] == 
                                                  #'group3_mixed_repressor')]


#data_oc_rna = data[(data['group'] == 'group_1_2_activator') | 
                     #(data['group'] == 'group3_repressor') | 
                     #(data['group'] == 'group3_mixed_repressor')]


# In[11]:


## Define: X_data (OC - open chromatin) & 
## y_data (gene expression)

def get_data(data):
    "input: dataframe with positive and negative correlations between oc and rna"
    "output: X_data, y_data for the model"

    #oc
    #start_column_x = 'X10.5_oc'
    #end_column_x = 'X21.0_oc'
    
    start_column_x = 'X10.5_oc1'   ## 60: X10.5_oc1
    end_column_x = 'X21.0_oc60'    ##  60: X21.0_oc60

    start_column_x = data.columns.get_loc(start_column_x)
    end_column_x = data.columns.get_loc(end_column_x)

    #rna
    start_column_y = 'X10.5_rna'      
    end_column_y = 'X21.0_rna'
    start_column_y = data.columns.get_loc(start_column_y)
    end_column_y = data.columns.get_loc(end_column_y)

    oc = data.iloc[:,start_column_x:end_column_x+1]
    rna = data.iloc[:,start_column_y:end_column_y+1]

    X_data = oc.values
    y_data = rna.values
    
    return (X_data, y_data)


# In[84]:


X_data_OC_, y_data_rna_ = get_data(data_oc_rna)


# In[85]:


X_data_OC_ = X_data_OC_.reshape(data_oc_rna.shape[0], num_cCREs_per_time_p, 8)


# In[86]:


### normalaized
X_data_OC = np.array(min_max_norm(X_data_OC_))
y_data_rna = np.array(min_max_norm(y_data_rna_))

### (rescaled) not normalized 
X_data_OC_not_norm = np.array(X_data_OC_)
y_data_rna_not_norm = np.array(y_data_rna_)


# In[87]:


def tensorxy(X_data_OC, y_data_rna):
    "input: X (oc) and y(rna) arrays "
    "output: X (oc) and y(rna) tensors"  
    
  
    seq_len = 8

    X_OC = X_data_OC.astype(float)



    seq_len_y = 8

    y_rna = y_data_rna.reshape(len(y_data_rna),1,seq_len_y)

    # convert into PyTorch tensors
    
    #X = torch.tensor(X, dtype=torch.float32)
    X_OC = torch.tensor(X_OC, dtype=torch.float32)

    #y = torch.tensor(y, dtype=torch.float32)
    y_rna = torch.tensor(y_rna, dtype=torch.float32)
    
   
    return(X_OC,y_rna)


# In[88]:


## Define X (Features) and y (target) as tensors

## normalized
X1,y1 = tensorxy(X_data_OC,y_data_rna)

## rescaled (not normalized)
X1_not_norm,y1_not_norm = tensorxy(X_data_OC_not_norm,
                                   y_data_rna_not_norm)


# ### NN model 

# In[89]:


batch_size_train = 4 ## can be changed to improve performance, larger--> runs faster. should be changed in train dataloader as well
batch_size_val = 1 
batch_size_test = 1 


# In[90]:


#split the data into training set (80%) and temporary set (20%)
def train_test_sets(X1, y1):
    
    X_train1, X_test1, y_train1, y_test1, train_indexes1, test_indexes1 = train_test_split(X1, y1, np.arange(len(X1)),test_size=0.2, random_state=42)

    #n_points = 105
    n_points = 8
    num_features = 60 ## 1 / 10 / 60

    X_train1 = X_train1.reshape(-1, num_features, n_points) 

    X_test1 = X_test1.reshape(-1, num_features, n_points) 

    #X_train1 = X_train1[:,0].reshape(X_train1.size()[0],num_features,n_points) # for OC feature only

    train_dataset1 = TensorDataset(X_train1[:,:,:n_points], y_train1[:,:,:n_points])

    train_loader1 = DataLoader(train_dataset1, batch_size=batch_size_train, shuffle=True)

    #X_test1 = X_test1[:,0].reshape(X_test1.size()[0],num_features,n_points) # use for OC feature only when num_features=1

    test_dataset1 = TensorDataset(X_test1[:,:,:n_points], y_test1[:,:,:n_points])

    test_loader1 = DataLoader(test_dataset1, batch_size=batch_size_test, shuffle=False)

    return train_loader1, test_loader1, X_test1, y_test1, test_indexes1 


# In[91]:


## normalized
train_loader1, test_loader1, X_test1_norm, y_test1_norm, test_indexes1_norm  = train_test_sets(X1, y1)

## rescaled (not normalized)
train_loader1_r, test_loader1_r, X_test1_r, y_test1_r, test_indexes1_r  = train_test_sets(X1_not_norm,y1_not_norm)


# In[92]:


torch.manual_seed(42)
num_samples = y1.shape[0]


time_steps = 8
num_features = 60    # 60 set number of features
batch_size = 4  # Define the batch size
loss_train_epoch = []


# Define network architacture to regress time series gene expression
#from open chromatin signals
class PolynomialRegression(nn.Module):
    
    def __init__(self, degree, num_features):
        super(PolynomialRegression, self).__init__()
        self.degree = degree
        
        self.poly = nn.Linear(num_features, degree, bias=True)
        self.rnn = nn.RNN(60, 30,bidirectional=True)   ### add / remove for linear  60,30
        self.relu = nn.ReLU(inplace=True) ## / remove for linear
        self.fc = nn.Linear(degree, 1) 
     
 
        
    def forward(self, x):
        
        batch_size, time_steps, num_features = x.size()
        output, hn  = self.rnn(x) ### add  / remove for linear
        x_poly = self.poly(output) ## add / remove for linear
        x_poly = self.relu(x_poly) ## add / remove for linear
        out = self.fc(x_poly)
        out = self.relu(out)  ## to keep the output positive / remove for linear
        
        return out

# Custom loss function to handle per-sample shape
class PerSampleMSELoss(nn.Module):
    def __init__(self):
        super(PerSampleMSELoss, self).__init__()
    

    def forward(self, input, target): 
        loss = torch.mean((input - target) ** 2, dim=1)  # Calculate the loss along the time_steps dimension
        return loss

# Initialize the model
degree = 10

warmup_steps = 1000  # Number of warm-up steps
initial_lr = 0.0001  #0.001  # Initial learning rate (warm-up learning rate)
desired_lr = 0.000001  #0.00001  # Desired learning rate


model = PolynomialRegression(degree, num_features).to(device)
model

# Define loss function and optimizer
criterion = PerSampleMSELoss()
optimizer = optim.Adam(model.parameters(), lr=initial_lr)

# Create a LambdaLR scheduler for warm-up
warmup_scheduler = LambdaLR(optimizer, lr_lambda=lambda step: step / warmup_steps if step < warmup_steps else 1)

# Training loop
num_epochs = 200    
for epoch in range(num_epochs):
 
    for i, (x_batch, y_batch) in enumerate(train_loader1):
        
        
        x_batch = x_batch.transpose(1, 2).requires_grad_().to(device)  
        y_batch = y_batch.transpose(1, 2).requires_grad_().to(device) 
        

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


# In[93]:


## run this cell before loading a Pytorch model

## save model ###  ## to save a Pytorch model use:
#torch.save(model, '')

degree = 10

class PolynomialRegression(nn.Module):
    def __init__(self, degree, num_features):
        super(PolynomialRegression, self).__init__()
        self.degree = degree
        
        
        self.poly = nn.Linear(num_features, degree, bias=True)
        self.rnn = nn.RNN(60, 30,bidirectional=True)   
        self.relu = nn.ReLU(inplace=True)  
        self.fc = nn.Linear(degree, 1)
     

    def forward(self, x):
        batch_size, time_steps, num_features = x.size()
        #x_poly = self.poly(x)  #add for linear NN
        output, hn  = self.rnn(x) 
        
        x_poly = self.poly(output)
        
        x_poly = self.relu(x_poly) 
        out = self.fc(x_poly)
        out = self.relu(out)  
        
        return out
    
# Custom loss function to handle per-sample shape

class PerSampleMSELoss(nn.Module):
    def __init__(self):
        super(PerSampleMSELoss, self).__init__()
    
        

    def forward(self, input, target):
        loss = torch.mean((input - target) ** 2, dim=1)  # Calculate the loss along the time_steps dimension
        
        return loss
criterion = PerSampleMSELoss()


# In[94]:


### save the first point of each gene from the rescaled table (NOT normalized)
rescaled_gene = y_test1_r
first_point_gene = [i[0][0] for i in rescaled_gene.tolist()]
#first_point_gene


# In[95]:


# Testing / inference
test_losses = []
all_pred_y = []
all_true_y = []

batch_size_test = 1
with torch.no_grad():
    
    for i, (x_batch, y_batch) in enumerate(test_loader1):
        
        x_batch = x_batch.transpose(1, 2).to(device)   
        y_batch = y_batch.transpose(1, 2).to(device)
        
        ## run this line if you run the training cell 
        #y_pred = model(x_batch)
        
        ## run this line to load the model 
        #model = torch.load('***.pth')

        
        model.eval()
        y_pred = model(x_batch)
        all_pred_y.append(y_pred)
        all_true_y.append(y_batch)

        # Compute the loss
       
        loss = criterion(y_pred, y_batch)
        print(f'Loss: {loss.mean().item():.4f}')
        test_losses.append(loss.mean().item())


# In[96]:


evaluated_G = rescaled_genes(all_pred_y, first_point_gene)


# In[97]:


# plot true gene expression versus predicted gene expression over time

num_points = 8
t =[10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.6, 21]

for i in range(0,len(all_true_y)):  
    plt.figure(figsize=(4,4))
    plt.scatter(t, all_true_y[i].cpu().numpy())
    plt.scatter(t, all_pred_y[i].cpu().numpy())
    
    corr_true_pred = corr(all_true_y[i].cpu().numpy(), rescaled_gene[i].reshape(1,1,8), [evaluated_G[i].cpu().numpy()])
    
    #plt.title(f'#{i} {data_oc_rna.iloc[test_indexes1_norm[i],[0]].values.tolist()[0],  data_oc_rna.iloc[test_indexes1_norm[i],[1]].values.tolist()[0]}')
    plt.title(f'#{i} {data_oc_rna.iloc[test_indexes1_norm[i],[0]].values.tolist()[0]}, \n corr true-pred: {corr_true_pred}')
    plt.xlabel('mouse PCDs', fontsize=18)
    
    plt.ylabel('G \n (gene expression)', fontsize=18)
    plt.ylabel('G', fontsize=18)
    
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim(10,22)
    plt.legend(['true G', 
                'predicted G'], 
               loc='upper center', bbox_to_anchor=(0.65, 1.5),fontsize=15 )


# ### Rescale back to gene expression data

# In[24]:


def rescaled_genes(all_pred_y, first_point_gene):
    
    evaluated_G = []
    
    for n, gene in enumerate(all_pred_y):
        
        first_point = first_point_gene[n]
        
        ### check if the first point is the min
        first_point_pred = gene[:,0]
        last_point_pred = gene[:,-1]
        
        if last_point_pred<first_point_pred:
            if first_point<1:
                first_point = 0
            else:
                first_point  = first_point - 1        
            
        evaluated_G.append(gene + first_point)   
        
    return evaluated_G


# In[98]:


corr_all = []


for i in range(0, len(y_test1_r)):
    plt.figure(figsize=(5,5))
    plt.scatter([10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.6, 21], rescaled_gene [i]) ## true gene expression from the rescaled table the test set
    plt.scatter([10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.6, 21], evaluated_G[i].cpu().numpy()) ## reconstructed
    
    corr_true_pred = corr(all_true_y[i].cpu().numpy(), rescaled_gene[i].reshape(1,1,8), [evaluated_G[i].cpu().numpy()])
    corr_all.append(corr_true_pred) 
    
    plt.title(f'#{i} {data_oc_rna.iloc[test_indexes1_r[i],[0]].values.tolist()[0]}, \n corr true-pred: {corr_true_pred}')
    plt.legend(['rescaled gene', 'evaluated gene'])
#filename = "brain_dev_many/60ccres/****.pdf"
#save_multi_image(filename)


# ## Correlation (Pearson's r) between the true and the predicted values

# In[25]:


## plot true-predicted correlation histogram

def corr(all_true_y, rescaled_gene, evaluated_G):

    corr_list = []
    true_list = []
    pred_list = []

    for i in range(0, len(all_true_y)):
        
        true_list.append(rescaled_gene[i].tolist()[0]) 
        pred_list.append([i[0] for i in evaluated_G[i][0].tolist()])  ## evaluated gene expression


    for i in range(0, len(pred_list)):
        r,p = sc.pearsonr(pred_list[i], true_list[i])
        corr_list.append(r)

    corr_list = [0 if math.isnan(x) else x for x in corr_list]
    
    
    return corr_list


# In[100]:


### computed correlations
print('mean correlation',  np.mean(corr_all))
print('std correlation',  np.std(corr_all))
print('median correlation',  np.median(corr_all))
print('IQR correlation',  iqr(corr_all))

plt.figure(figsize=(6.3,6.3))
plt.hist([i[0] for i in corr_all], histtype = 'stepfilled', color='royalblue')
plt.xticks(fontsize=13)
plt.xlabel('true - predicted G \n correlation', fontsize=13)
plt.yticks(fontsize=13)
plt.ylabel('counts', fontsize=13)
plt.xticks([i for i in range(-1,2,1)])
corr_table = pd.DataFrame(corr_all)
#corr_table.to_csv('***.csv')
#filename = "brain_dev_many/60ccres/****.pdf"
#save_multi_image(filename)


# In[101]:


ax = plt.figure(figsize = (6,6))
true_list = []
pred_list = []

for i in range(0, len(y_test1_r)):
    true_list.append(rescaled_gene[i].tolist()[0])  
    pred_list.append([i[0] for i in evaluated_G[i][0].tolist()])

for i in range(0, len(evaluated_G)):

    plt.scatter(true_list[i], pred_list[i], c = 'blue', s = 1, alpha=0.3)
    plt.plot([i for i in range(0, 12)], [i for i in range(0, 12)], '--', c = 'black')
    plt.xticks(fontsize = 13)
    plt.yticks(fontsize = 13)
    plt.xlabel('True values', fontsize = 15)
    plt.ylabel('Predicted values', fontsize = 15)

# Setting the axis range for both axes
plt.xticks([i for i in range(0,13,2)])
plt.yticks([i for i in range(0,13,2)])

lt = []
for i in true_list:
    lt.extend(i)
lp = []
for i in pred_list:
    lp.extend(i)

d = {'true': lt, 'pred': lp}
my_table = pd.DataFrame(d)

#my_table.to_csv('****.csv')
#filename = "****.pdf"
#save_multi_image(filename)


# ## MSE

# In[102]:


## calculate mean and std of MSE
plt.figure(figsize=(2,2))

MSE_list = []

MSE_list = [np.sum(np.power(np.absolute(i),2))/8 for n,i in enumerate(np.subtract(true_list,pred_list))]

counts, bins = np.histogram(MSE_list , bins=50)  # Calculate counts and bin edges
relative_freq = counts / counts.sum()  # Calculate relative frequency

print(np.median(np.array(MSE_list)))
print(np.std(np.array(MSE_list)))


plt.bar(bins[:-1], relative_freq, width=np.diff(bins), alpha=0.7)   ## plot the MSE with fraction in the y axis
plt.ylim(0,1)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim(0,6)
plt.xlabel('MSE', fontsize=18)
plt.ylabel('Fraction', fontsize=18)    ### was counts

d = {'MSE': MSE_list}
my_table = pd.DataFrame(d)
#my_table.to_csv('****.csv')
#filename = "****.pdf"
#save_multi_image(filename)


# In[92]:


plt.figure(figsize=(2,2))
plt.hist(my_MSE['MAE'].tolist(), bins=50)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(0,10)
plt.xlabel('MSE', fontsize=16)
plt.ylabel('counts', fontsize=16)

