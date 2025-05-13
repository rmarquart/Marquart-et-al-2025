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
iceFloeCategoryAverages = np.array(df['iceFloeCategoryAverages'])
print(f"iceFloeCategoryAverages = {iceFloeCategoryAverages}")

# Directory where files are stored
directory = f'/scratch/mrqrut001/simulation/postdoc/Python_real/Previous index values/2022/dt=450s/h1/T{num}'

# Ensure directory path
os.chdir(directory)

# Identify positions of non-NaN values in iceFloeCategoryAverages
non_nan_positions = [i for i, val in enumerate(iceFloeCategoryAverages) if not np.isnan(val)]
print(f"non_nan_positions = {non_nan_positions}")

# Check if all required files exist
all_files_exist_h1 = all(
    os.path.exists(os.path.join(directory, f'prev_index_values_{pos}.pkl'))
    for pos in non_nan_positions
)

print(f"all_files_exist_h1 = {all_files_exist_h1}")

# If all files exist, skip the loop
if all_files_exist_h1:
    print("All necessary files exist. Skipping the loop.")
else:
    print("Some files are redundant or missing, processing loop...")
    for i, category in enumerate(iceFloeCategoryAverages):
        filename_to_remove = f'prev_index_values_{i}.pkl'

        # Initialize `loaded_variables` to avoid undefined access
        loaded_variables = None

        # If category is NaN, remove existing file
        if np.isnan(category):
            if os.path.exists(filename_to_remove):
                print(f"Removing {filename_to_remove}")
                os.remove(filename_to_remove)

        # Otherwise, handle non-NaN categories
        else:
            filename = f'prev_index_values_{i}.pkl'
            alternative_filename1 = f'prev_index_values_{i+1}.pkl'
            alternative_filename2 = f'prev_index_values_{i-1}.pkl'
            alternative_filename3 = f'prev_index_values_{i+2}.pkl'
            alternative_filename4 = f'prev_index_values_{i-2}.pkl'

            # Check if filename exists; otherwise, use alternatives
            if not os.path.exists(filename):
                if os.path.exists(alternative_filename1):
                    shutil.copy(alternative_filename1, filename)
                    print(f"Copied {alternative_filename1} to {filename}")

                    with open(alternative_filename1, 'rb') as f:
                        loaded_variables = pickle.load(f)

                    with open(filename, 'wb') as f:
                        pickle.dump(loaded_variables, f)

                    print(f"Saving file: {filename}")
                    print(f"'index' variable for {filename}: {loaded_variables['index']}")

                elif os.path.exists(alternative_filename2):
                    shutil.copy(alternative_filename2, filename)
                    print(f"Copied {alternative_filename2} to {filename}")

                    with open(alternative_filename2, 'rb') as f:
                        loaded_variables = pickle.load(f)

                    with open(filename, 'wb') as f:
                        pickle.dump(loaded_variables, f)

                    print(f"Saving file: {filename}")
                    print(f"'index' variable for {filename}: {loaded_variables['index']}")

                elif os.path.exists(alternative_filename3):
                    shutil.copy(alternative_filename3, filename)
                    print(f"Copied {alternative_filename3} to {filename}")

                    with open(alternative_filename3, 'rb') as f:
                        loaded_variables = pickle.load(f)

                    with open(filename, 'wb') as f:
                        pickle.dump(loaded_variables, f)

                    print(f"Saving file: {filename}")
                    print(f"'index' variable for {filename}: {loaded_variables['index']}")

                elif os.path.exists(alternative_filename4):
                    shutil.copy(alternative_filename4, filename)
                    print(f"Copied {alternative_filename4} to {filename}")

                    with open(alternative_filename4, 'rb') as f:
                        loaded_variables = pickle.load(f)

                    with open(filename, 'wb') as f:
                        pickle.dump(loaded_variables, f)

                    print(f"Saving file: {filename}")
                    print(f"'index' variable for {filename}: {loaded_variables['index']}")

                else:
                    print(f"Neither {alternative_filename1}, {alternative_filename2}, {alternative_filename3} nor {alternative_filename4} exists.")

# Define the path to the existing 'storage.csv' file
csv_filepath = "/scratch/mrqrut001/simulation/postdoc/Python_real/Storage csv/T{}/storage.csv".format(num)

# Read the existing CSV file into a DataFrame
df = pd.read_csv(csv_filepath)

# Define the column name and value
new_column = 'all_files_exist_h1'
new_value = all_files_exist_h1

# Check if the column exists
df[new_column] = new_value

# Save the updated DataFrame back to the CSV file
df.to_csv(csv_filepath, index=False)
