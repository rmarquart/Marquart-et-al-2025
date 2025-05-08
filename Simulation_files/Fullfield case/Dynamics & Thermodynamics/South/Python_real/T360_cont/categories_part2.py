# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 15:13:35 2023

@author: Ma-re Admin
"""

import numpy as np
import os
import pandas as pd
import math

# Change directory to 'Storage h' folder
os.chdir('/scratch/mrqrut001/simulation/postdoc/Python_real/')

# Load the exported variable
numCategories_h1 = int(os.environ.get("numCategories_h1"))
numCategories_h2 = int(os.environ.get("numCategories_h2"))
num = str(os.environ.get("num"))

# Define the path to the 'storage.csv' file
csv_filepath = "/scratch/mrqrut001/simulation/postdoc/Python_real/Storage csv/T{}/storage.csv".format(num)

# Read the CSV file into a DataFrame
df = pd.read_csv(csv_filepath)

# Access the values from the DataFrame
iceFloeCategoryAverages = np.array(df['iceFloeCategoryAverages'])
greaseIceCategoryAverages = np.array(df['greaseIceCategoryAverages'])
numCellsIceFloe = np.array(df['numCellsIceFloe'])
numCellsGreaseIce = np.array(df['numCellsGreaseIce'])
all_files_exist_h1 = df['all_files_exist_h1'].iloc[0]
all_files_exist_h2 = df['all_files_exist_h2'].iloc[0]

# if all_files_exist_h1:
print(f"all_files_exist_h1 = {all_files_exist_h1}")
# Initialize an empty list to store hi_cat values
h1_cat = []

for i, category in enumerate(iceFloeCategoryAverages):
    if np.isnan(category):
        h1_cat.append(np.nan)  # Optionally, add NaN directly to h1_cat
    else:
        column_name = f'h1_cat{i}'  # Ensure category is converted to int if needed
        if column_name in df.columns:
            h1_cat_values = df[column_name].iloc[0]
            h1_cat.append(h1_cat_values)
        else:
            h1_cat.append(np.nan)

# Convert h1_cat to a NumPy array if needed
h1_cat = np.array(h1_cat)

# Define different offsets for each category
offsetsIceFloe = np.array([h1_cat[0]-iceFloeCategoryAverages[0], h1_cat[1]-iceFloeCategoryAverages[1], h1_cat[2]-iceFloeCategoryAverages[2], h1_cat[3]-iceFloeCategoryAverages[3], h1_cat[4]-iceFloeCategoryAverages[4], h1_cat[5]-iceFloeCategoryAverages[5]])

# Ensure that offsetsIceFloe and offsetsGreaseIce have the same length
max_length = len(df)
offsetsIceFloe = np.pad(offsetsIceFloe, (0, max_length - len(offsetsIceFloe)), mode='constant', constant_values=np.nan)

numCellsIceFloeTotal = numCellsIceFloe[0]+numCellsIceFloe[1]+numCellsIceFloe[2]+numCellsIceFloe[3]+numCellsIceFloe[4]+numCellsIceFloe[5]

dhIceFloe_cat0 = offsetsIceFloe[0]
dhIceFloe_cat1 = offsetsIceFloe[1]
dhIceFloe_cat2 = offsetsIceFloe[2]
dhIceFloe_cat3 = offsetsIceFloe[3]
dhIceFloe_cat4 = offsetsIceFloe[4]
dhIceFloe_cat5 = offsetsIceFloe[5]

# Replace NaN values with zero
offsetsIceFloe[np.isnan(offsetsIceFloe)] = 0

for i, category in enumerate(iceFloeCategoryAverages):
    print(f'Offsets for ice floe category {i}:')
    print(f'   Offset: {offsetsIceFloe[i]:.6f}')

del h1_cat

# if all_files_exist_h2:
print(f"all_files_exist_h2 = {all_files_exist_h2}")

# Initialize an empty list to store hi_cat values
h2_cat = []

for i, category in enumerate(greaseIceCategoryAverages):
    if np.isnan(category):
        h2_cat.append(np.nan)  # Optionally, add NaN directly to h1_cat
    else:
        column_name = f'h2_cat{i}'  # Ensure category is converted to int if needed
        if column_name in df.columns:
            h2_cat_values = df[column_name].iloc[0]
            h2_cat.append(h2_cat_values)
        else:
            h2_cat.append(np.nan)

# Convert h2_cat to a NumPy array if needed
h2_cat = np.array(h2_cat)

# Define different offsets for each category
offsetsGreaseIce = np.array([h2_cat[0]-greaseIceCategoryAverages[0], h2_cat[1]-greaseIceCategoryAverages[1], h2_cat[2]-greaseIceCategoryAverages[2], h2_cat[3]-greaseIceCategoryAverages[3], h2_cat[4]-greaseIceCategoryAverages[4], h2_cat[5]-greaseIceCategoryAverages[5]])

# Ensure that offsetsIceFloe and offsetsGreaseIce have the same length
max_length = len(df)
offsetsGreaseIce = np.pad(offsetsGreaseIce, (0, max_length - len(offsetsGreaseIce)), mode='constant', constant_values=np.nan)

dhGreaseIce_cat0 = offsetsGreaseIce[0]
dhGreaseIce_cat1 = offsetsGreaseIce[1]
dhGreaseIce_cat2 = offsetsGreaseIce[2]
dhGreaseIce_cat3 = offsetsGreaseIce[3]
dhGreaseIce_cat4 = offsetsGreaseIce[4]
dhGreaseIce_cat5 = offsetsGreaseIce[5]

# Replace NaN values with zero
offsetsGreaseIce[np.isnan(offsetsGreaseIce)] = 0

for i, category in enumerate(greaseIceCategoryAverages):
     print(f'Offsets for grease ice category {i}:')
     print(f'   Offset: {offsetsGreaseIce[i]:.6f}')

del h2_cat

# Define the path to the existing 'storage.csv' file
csv_filepath = "/scratch/mrqrut001/simulation/postdoc/Python_real/Storage csv/T{}/storage.csv".format(num)

# Read the existing CSV file into a DataFrame
df = pd.read_csv(csv_filepath)

value = 'offsetsIceFloe'
if value in df.iloc[0]:
    # Update the columns with new values
    df[value] = offsetsIceFloe
else:
    # You can add the new variables as new columns in the DataFrame
    df = pd.concat([df,pd.DataFrame({value : offsetsIceFloe})],axis=1)

value = 'offsetsGreaseIce'
if value in df.iloc[0]:
    # Update the columns with new values
    df[value] = offsetsGreaseIce
else:
    # You can add the new variables as new columns in the DataFrame
    df = pd.concat([df,pd.DataFrame({value : offsetsGreaseIce})],axis=1)

# Save the updated DataFrame back to the same CSV file
df.to_csv(csv_filepath, index=False)