# Author: Eve Wattenberg
# Random Forest regression

# Imports
import numpy as np
import pandas as pd
import argparse
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score


# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputfile", type=str,
                help="Input file with paired RNA and OC values")
parser.add_argument("-p", "--predictfile", type=str,
                help="File to write predictions to")  
parser.add_argument("-t", "--timepoints", type=str, default=105,
                help="Number of time points in the data")  
parser.add_argument("-m", "--modalities", type=str, default=1,
                help="Number of input data modalities")  
parser.add_argument("-s", "--split", type=str, default=0.2,
                help="Fraction of input to use as test set")                
parser.add_argument("-r", "--randomseed", type=str, default=1941,
                help="Random seed")        
args = parser.parse_args()

# Set random seed
np.random.seed(args.randomseed)

# Read in input/output file.
raw = pd.read_csv(args.inputfile, sep="\t")
raw['combined_index'] = raw.apply(lambda row: f"{row['gene_id']}-{row['region']}-{row['cCRE_id']}", axis=1)
raw.index = raw["combined_index"]

# Select RNA cols
rna_start = (args.timepoints * args.modalities) + 3
rna_end = rna_start + args.timepoints
all_rna = raw.iloc[:, rna_start:rna_end]
all_rna.index = raw.index

# Select input cols
all_input = raw.iloc[:, 3:rna_start]
all_input.index = raw.index

# Split into training and test sets
input_train, input_test, rna_train, rna_test = train_test_split(all_input, all_rna, test_size=args.split, random_state=42) # 30% data as test set


# Model!
rf_model = RandomForestRegressor(n_estimators=100) # arbitrary number of trees

# Train the model
rf_model.fit(input_train, rna_train)

# Test the model
predictions = rf_model.predict(input_test)

# Compute errors
mae = mean_absolute_error(rna_test, predictions)
mse = mean_squared_error(rna_test, predictions)
r2 = r2_score(rna_test, predictions)
print("Mean Absolute Error:", mae)
print("Mean Squared Error:", mse)
print("R^2 Score:", r2)

# Save predictions
predictions = pd.DataFrame(predictions)
predictions.index = rna_test.index
predictions.to_csv(args.predictfile, sep="\t")



