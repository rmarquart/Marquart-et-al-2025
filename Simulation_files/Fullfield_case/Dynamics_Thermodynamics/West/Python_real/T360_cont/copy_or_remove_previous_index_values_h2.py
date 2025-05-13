# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 14:57:11 2023

@author: Ma-re Admin
"""
import pickle 
import os
import pandas as pd
import numpy as np
import shutil

num = str(os.environ.get("num"))

# Define the path to the 'storage.csv' file
csv_filepath = "/scratch/mrqrut001/simulation/postdoc/Python_real/Storage csv/T{}/storage.csv".format(num)

# Read the CSV file into a DataFrame
df = pd.read_csv(csv_filepath)

# Access the values from the DataFrame
greaseIceCategoryAverages = np.array(df['greaseIceCategoryAverages'])
print(f"greaseIceCategoryAverages = {greaseIceCategoryAverages}")

# Directory where files are stored
directory = f'/scratch/mrqrut001/simulation/postdoc/Python_real/Previous index values/2022/dt=450s/h2/T{num}'

# Ensure directory path
os.chdir(directory)

# Identify positions of NaN values in greaseIceCategoryAverages
nan_positions = [i for i, val in enumerate(greaseIceCategoryAverages) if np.isnan(val)]
print(f"NaN positions = {nan_positions}")

# **Step 1: Remove all .pkl files corresponding to NaN positions**
for i in nan_positions:
    filename_to_remove = os.path.join(directory, f'prev_index_values_{i}.pkl')
    if os.path.exists(filename_to_remove):
        print(f"Removing {filename_to_remove}")
        os.remove(filename_to_remove)  # Ensure the file gets removed

# **Step 2: Identify non-NaN positions and check if all required files exist**
non_nan_positions = [i for i in range(len(greaseIceCategoryAverages)) if i not in nan_positions]
all_files_exist_h2 = all(
    os.path.exists(os.path.join(directory, f'prev_index_values_{pos}.pkl'))
    for pos in non_nan_positions
)

print(f"all_files_exist_h2 = {all_files_exist_h2}")

# **Step 3: If required files are missing, copy them from adjacent ones**
if not all_files_exist_h2:
    print("Some required files are missing, processing replacements...")

    for i in non_nan_positions:
        filename = os.path.join(directory, f'prev_index_values_{i}.pkl')

        # If the file does not exist, attempt to copy from adjacent files
        if not os.path.exists(filename):
            for offset in [1, -1, 2, -2]:  # Prioritize closest files
                alternative_filename = os.path.join(directory, f'prev_index_values_{i + offset}.pkl')
                
                if os.path.exists(alternative_filename):
                    shutil.copy(alternative_filename, filename)
                    print(f"Copied {alternative_filename} to {filename}")

                    with open(alternative_filename, 'rb') as f:
                        loaded_variables = pickle.load(f)

                    with open(filename, 'wb') as f:
                        pickle.dump(loaded_variables, f)

                    print(f"Saved file: {filename}")
                    print(f"'index' variable for {filename}: {loaded_variables['index']}")
                    break
            else:
                print(f"No suitable adjacent file found for {filename}")

# Define the path to the existing 'storage.csv' file
csv_filepath = "/scratch/mrqrut001/simulation/postdoc/Python_real/Storage csv/T{}/storage.csv".format(num)

# Read the existing CSV file into a DataFrame
df = pd.read_csv(csv_filepath)

# Define the column name and value
new_column = 'all_files_exist_h2'
new_value = all_files_exist_h2

# Check if the column exists
df[new_column] = new_value

# Save the updated DataFrame back to the CSV file
df.to_csv(csv_filepath, index=False)
