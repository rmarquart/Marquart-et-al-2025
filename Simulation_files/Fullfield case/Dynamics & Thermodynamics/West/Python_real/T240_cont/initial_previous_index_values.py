# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 14:57:11 2023

@author: Ma-re Admin
"""
import pickle 
import os
import pandas as pd
import numpy as np

#####################################  h1  #####################################

# Load the exported variable
index = int(os.environ.get("index"))
num = str(os.environ.get("num"))

# Define the path to the 'storage.csv' file
csv_filepath = "/scratch/mrqrut001/simulation/postdoc/Python_real/Storage csv/T{}/storage.csv".format(num)

# Read the CSV file into a DataFrame
df = pd.read_csv(csv_filepath)

# Access the values from the DataFrame
iceFloeCategoryAverages = np.array(df['iceFloeCategoryAverages'])
print(f"iceFloeCategoryAverages = {iceFloeCategoryAverages}")

# Change directory to 'Saved variables' folder
os.chdir('/scratch/mrqrut001/simulation/postdoc/Python_real/Saved variables/2022/dt=450s/h1/fullfield')

# Define the values for i and j
i_values = [14, 16, 18]  # Add more values to i_values if needed
j_values = range(len(i_values))  # Corresponding j_values for categories

# Initialize an empty dictionary to store the loaded variables
loaded_variables = {}

# Track the count of non-NaN elements found
nonzero_count = 0

for idx in range(0, len(iceFloeCategoryAverages)):
    if np.isnan(iceFloeCategoryAverages[idx]):
        print(f"Category {idx} is empty, thus no need to for a previous time step.")
    else:
        # Increment the non-NaN count
        nonzero_count += 1

        # Check if there is a corresponding i_value for the current nonzero_count
        if nonzero_count <= len(i_values):
            # Select the nth value in i_values for the nth nonzero element
            i = i_values[nonzero_count - 1]  # 0-based index
            i_str = f'{i:02d}'

            # Generate the modified filename
            filename = f'saved_variables_0.{i_str}m_new.pkl'

            # Open the file with the modified filename
            with open(filename, 'rb') as f:
                # Name the loaded variable based on i
                loaded_var_name = f'loaded_variables_h{i_str}'
                loaded_variables[loaded_var_name] = pickle.load(f)

            # Check if the loaded variable has a 'hi' key
            if 'hi' in loaded_variables[loaded_var_name]:
                # Create a new variable like h01 for loaded_variables_h01['hi']
                h_variable_name = f'h{i_str}'
                globals()[h_variable_name] = loaded_variables[loaded_var_name]['hi']

# Collect initial thickness values if they exist
initial_thicknesses_h1 = []
for i in i_values:
    h_variable_name = f'h{i:02d}'
    if h_variable_name in globals():
        initial_thicknesses_h1.append(globals()[h_variable_name][index])

# Print the initial thicknesses with the specified format
for j, thickness in zip(j_values, initial_thicknesses_h1):
    print(f'Category {j + 1}: initial thickness h1 = {thickness:.6f}m')
            
values_at_index = {}

prev_index = int(index) - 1

# Define a function to retrieve values at a specific index in arrays within a dictionary
def get_values_at_index(dictionary, prev_index):
    values_at_index = {key: array[prev_index] for key, array in dictionary.items() if len(array) > prev_index}
    values_at_index["index"] = prev_index
    return values_at_index

loaded_variables_list = [loaded_variables[f'loaded_variables_h{i:02d}'] for i in i_values]

# Change directory to 'Previous index values' folder
os.chdir(f'/scratch/mrqrut001/simulation/postdoc/Python_real/Previous index values/2022/dt=450s/h1/T{num}')

for j, loaded_variables in enumerate(loaded_variables_list):
    values_at_index = get_values_at_index(loaded_variables, prev_index)
    filename = f'prev_index_values_{j + 1}.pkl'
    print(f"Saving file: {filename}")
    with open(filename, 'wb') as f:
        pickle.dump(values_at_index, f)

# Define the path for the new CSV file
csv_filepath = "/scratch/mrqrut001/simulation/postdoc/Python_real/Storage csv/T{}/storage.csv".format(num)

# Read the existing CSV file into a DataFrame
df = pd.read_csv(csv_filepath)

# You can add the new variables as new columns in the DataFrame
df = pd.concat([df,pd.DataFrame({'initial_thicknesses_h1': initial_thicknesses_h1})],axis=1)

# Save the updated DataFrame back to the same CSV file
df.to_csv(csv_filepath, index=False)

print("")

# #####################################  h2  #####################################

# Load the exported variable
index = int(os.environ.get("index"))
num = str(os.environ.get("num"))

# Define the path to the 'storage.csv' file
csv_filepath = "/scratch/mrqrut001/simulation/postdoc/Python_real/Storage csv/T{}/storage.csv".format(num)

# Read the CSV file into a DataFrame
df = pd.read_csv(csv_filepath)

# Access the values from the DataFrame
greaseIceCategoryAverages = np.array(df['greaseIceCategoryAverages'])
print(f"greaseIceCategoryAverages = {greaseIceCategoryAverages}")

# Change directory to 'Saved variables' folder
os.chdir('/scratch/mrqrut001/simulation/postdoc/Python_real/Saved variables/2022/dt=450s/h2/fullfield')

# Define the values for i and j
i_values = [0]  # Add more values to i_values if needed
j_values = range(len(i_values))  # Corresponding j_values for categories

# Initialize an empty dictionary to store the loaded variables
loaded_variables = {}

# Track the count of non-NaN elements found
nonzero_count = 0

for idx in range(0, len(greaseIceCategoryAverages)):
    if np.isnan(greaseIceCategoryAverages[idx]):
        print(f"Category {idx} is empty, thus no need to for a previous time step.")
    else:
        # Increment the non-NaN count
        nonzero_count += 1

        # Check if there is a corresponding i_value for the current nonzero_count
        if nonzero_count <= len(i_values):
            # Select the nth value in i_values for the nth nonzero element
            i = i_values[nonzero_count - 1]  # 0-based index
            i_str = f'{i:02d}'

            # Generate the modified filename
            filename = f'saved_variables_{i_str}_grease_ice.pkl'

            # Open the file with the modified filename
            with open(filename, 'rb') as f:
                # Name the loaded variable based on i
                loaded_var_name = f'loaded_variables_h{i_str}'
                loaded_variables[loaded_var_name] = pickle.load(f)

            # Check if the loaded variable has a 'hi' key
            if 'hi' in loaded_variables[loaded_var_name]:
                # Create a new variable like h01 for loaded_variables_h01['hi']
                h_variable_name = f'h{i_str}'
                globals()[h_variable_name] = loaded_variables[loaded_var_name]['hi']

# Collect initial thickness values if they exist
initial_thicknesses_h2 = []
for i in i_values:
    h_variable_name = f'h{i:02d}'
    if h_variable_name in globals():
        initial_thicknesses_h2.append(globals()[h_variable_name][index])

# Print the initial thicknesses with the specified format
for j, thickness in zip(j_values, initial_thicknesses_h2):
    print(f'Category {j + 1}: initial thickness h2 = {thickness:.6f}m')
            
values_at_index = {}

prev_index = int(index) - 1

# Define a function to retrieve values at a specific index in arrays within a dictionary
def get_values_at_index(dictionary, prev_index):
    values_at_index = {key: array[prev_index] for key, array in dictionary.items() if len(array) > prev_index}
    values_at_index["index"] = prev_index
    return values_at_index

loaded_variables_list = [loaded_variables[f'loaded_variables_h{i:02d}'] for i in i_values]

# Change directory to 'Previous index values' folder
os.chdir(f'/scratch/mrqrut001/simulation/postdoc/Python_real/Previous index values/2022/dt=450s/h2/T{num}')

for j, loaded_variables in enumerate(loaded_variables_list):
    values_at_index = get_values_at_index(loaded_variables, prev_index)
    filename = f'prev_index_values_{j + 1}.pkl'
    print(f"Saving file: {filename}")
    with open(filename, 'wb') as f:
        pickle.dump(values_at_index, f)

# Define the path for the new CSV file
csv_filepath = "/scratch/mrqrut001/simulation/postdoc/Python_real/Storage csv/T{}/storage.csv".format(num)

# Read the existing CSV file into a DataFrame
df = pd.read_csv(csv_filepath)

# You can add the new variables as new columns in the DataFrame
df = pd.concat([df,pd.DataFrame({'initial_thicknesses_h2': initial_thicknesses_h2})],axis=1)

# Save the updated DataFrame back to the same CSV file
df.to_csv(csv_filepath, index=False)