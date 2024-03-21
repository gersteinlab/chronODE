# Author: ESW
# Random Forest regression; REFORMATTED after biorXiv submission

# Imports
import numpy as np
import pandas as pd
import sys
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

# Dear Past Eve: WHAT IS GOING ON WITH THE CLI INPUTS!!!

np.random.seed(1941)
# Read in input/output file.
raw = pd.read_csv(sys.argv[1], sep="\t")
raw['combined_index'] = raw.apply(lambda row: f"{row['gene_id']}-{row['region']}-{row['cCRE_id']}", axis=1)
raw.index = raw["combined_index"]
# Select RNA cols
all_rna = raw.iloc[:, 108:213]
all_rna.index = raw.index
print("All rna shape: ", all_rna.shape)
print(all_rna.head)

# Select input cols: oc, k27ac, and CORRELATION   
# all_input = raw.drop(columns=['ccre_id', 'gene_id', 'cCRE_id'])
# all_input = pd.DataFrame(all_input.filter(regex='^(?!.*rna).*$').values)
all_input = raw.iloc[:, 3:108]
all_input.index = raw.index
#all_input = all_input.drop([0, 1], axis=1)
print("All input shape: ", all_input.shape)
print(all_input.head)

# Split into training and train_test_split
input_train, input_test, rna_train, rna_test = train_test_split(all_input, all_rna, test_size=0.2, random_state=42) # 30% data as test set


# Model!
rf_model = RandomForestRegressor(n_estimators=100) # arbitrary number of trees

# Train
rf_model.fit(input_train, rna_train)
print("Training shape: ", rna_train.shape)
print(rna_train.head)

# Test
predictions = rf_model.predict(input_test)
print("Predictions shape:", predictions.shape)
print(predictions[:2])

mae = mean_absolute_error(rna_test, predictions)
mse = mean_squared_error(rna_test, predictions)
r2 = r2_score(rna_test, predictions)

print("Mean Absolute Error:", mae)
print("Mean Squared Error:", mse)
print("R^2 Score:", r2)

# Save predictions
predictions = pd.DataFrame(predictions)
predictions.index = rna_test.index
predictions.to_csv(sys.argv[2], sep="\t")

input_train.to_csv(sys.argv[3], sep="\t")
rna_train.to_csv(sys.argv[4], sep="\t")