# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 15:06:58 2023

@author: Ma-re Admin
"""

import pandas as pd
import numpy as np
import os

# Load the exported variable
num = str(os.environ.get("num"))

# Change directory to 'Storage h' folder
os.chdir("/scratch/mrqrut001/simulation/postdoc/Python_real/Storage eta/T{}".format(num))

# Load the 'eta' file
filename = f'eta{num}.npy'
eta = np.load(filename)

# Change directory to 'Storage h' folder
os.chdir("/scratch/mrqrut001/simulation/postdoc/Python_real/Storage h/T{}".format(num))

# Load the 'h' file
filename = f'h{num}.npy'
h = np.load(filename)

# Change directory to 'Storage alpha' folder
os.chdir("/scratch/mrqrut001/simulation/postdoc/Python_real/Storage alpha/T{}".format(num))

# Load the 'alpha' file
filename = f'alpha{num}.npy'
alpha = np.load(filename)

# Load the exported variable
numCategories_h1 = int(os.environ.get("numCategories_h1"))
numCategories_h2 = int(os.environ.get("numCategories_h2"))

# Initialize lists to hold thickness values and original positions for each category
iceFloeCategories = [[] for _ in range(numCategories_h1)]
iceFloePositions = [[] for _ in range(numCategories_h1)]
greaseIceCategories = [[] for _ in range(numCategories_h2)]
greaseIcePositions = [[] for _ in range(numCategories_h2)]

# Initialize arrays to hold the number of cells and averages for each category
numCellsIceFloe = np.zeros(numCategories_h1, dtype=int)
numCellsGreaseIce = np.zeros(numCategories_h2, dtype=int)
iceFloeCategoryAverages = np.zeros(numCategories_h1)
greaseIceCategoryAverages = np.zeros(numCategories_h2)

# Initialize lists to hold the boundaries of each category for ice floes and grease ice
iceFloeCategoryBounds = [[] for _ in range(numCategories_h1)]
greaseIceCategoryBounds = [[] for _ in range(numCategories_h2)]

# Define logical indices for ice floe and grease ice
iceFloeAlphaIndices = (alpha >= 0.8) & (eta > 1e4)
greaseIceAlphaIndices = ~iceFloeAlphaIndices

# Define fixed thickness categories with a small epsilon to avoid overlap for ice floes
epsilon = np.finfo(float).eps
iceFloe_category_bounds = [
    (0, 0.10 - epsilon),
    (0.10, 0.12 - epsilon),
    (0.12, 0.14 - epsilon),
    (0.14, 0.16 - epsilon),
    (0.16, 0.18 - epsilon),
    (0.18, 0.28)
]

# Define fixed thickness categories for grease ice
greaseIce_category_bounds = [
    (0, 0.06 - epsilon),
    (0.06, 0.09 - epsilon),
    (0.09, 0.12 - epsilon),
    (0.12, 0.15 - epsilon),
    (0.15, 0.18 - epsilon),
    (0.18, 0.24)
]

# Loop through each category for ice floes and store bounds and statistics
for i, (lower_bound, upper_bound) in enumerate(iceFloe_category_bounds):
    # Find indices within the current category bounds
    iceFloeCategoryIndices = np.where(iceFloeAlphaIndices & (h >= lower_bound) & (h < upper_bound))
    
    # Store indices and calculate number of cells in this category
    iceFloePositions[i] = iceFloeCategoryIndices[0]
    numCellsIceFloe[i] = len(iceFloePositions[i])
    
    # Print the number of ice floe category elements
    print(f'Ice Floe Category {i} ({round(lower_bound, 2)}-{round(upper_bound, 2)}): number of cells = {numCellsIceFloe[i]}')
    
    # Calculate the average thickness for the category, if not empty
    if numCellsIceFloe[i] == 0:
        iceFloeCategoryAverages[i] = np.nan
    else:
        iceFloeCategoryAverages[i] = np.mean(h[iceFloeCategoryIndices])
    
    # Store the category bounds
    iceFloeCategoryBounds[i] = [lower_bound, upper_bound]

# Display the categorized data for ice floes
for i in range(numCategories_h1):
    print(f'Ice Floe Category {i}:')
    print(f'   Average Thickness: {iceFloeCategoryAverages[i]:.6f}')

print()

# Loop through each category for grease ice and store bounds and statistics
for i, (lower_bound, upper_bound) in enumerate(greaseIce_category_bounds):
    # Find indices within the current category bounds
    greaseIceCategoryIndices = np.where(greaseIceAlphaIndices & (h >= lower_bound) & (h < upper_bound))
    
    # Store indices and calculate number of cells in this category
    greaseIcePositions[i] = greaseIceCategoryIndices[0]
    numCellsGreaseIce[i] = len(greaseIcePositions[i])
    
    # Print the number of grease ice category elements
    print(f'Grease Ice Category {i} ({round(lower_bound, 3)}-{round(upper_bound, 3)}): number of cells = {numCellsGreaseIce[i]}')
    
    # Calculate the average thickness for the category, if not empty
    if numCellsGreaseIce[i] == 0:
        greaseIceCategoryAverages[i] = np.nan
    else:
        greaseIceCategoryAverages[i] = np.mean(h[greaseIceCategoryIndices])
    
    # Store the category bounds
    greaseIceCategoryBounds[i] = [lower_bound, upper_bound]

# Display the categorized data for grease ice
for i in range(numCategories_h2):
    print(f'Grease Ice Category {i}:')
    print(f'   Average Thickness: {greaseIceCategoryAverages[i]:.6f}')

test = "test"

df = pd.DataFrame({'test': [test]})

# Define the path for the new CSV file
csv_filepath = "/scratch/mrqrut001/simulation/postdoc/Python_real/Storage csv/T{}/storage.csv".format(num)

# Save the updated DataFrame back to the same CSV file
df.to_csv(csv_filepath, index=False)

# Read the existing CSV file into a DataFrame
df = pd.read_csv(csv_filepath)

# You can add the new variables as new columns in the DataFrame
df = pd.concat([df,pd.DataFrame({'iceFloeCategoryAverages': iceFloeCategoryAverages}),pd.DataFrame({'greaseIceCategoryAverages': greaseIceCategoryAverages}),pd.DataFrame({'numCellsIceFloe': numCellsIceFloe}),pd.DataFrame({'numCellsGreaseIce': numCellsGreaseIce})],axis=1)

# Save the updated DataFrame back to the same CSV file
df.to_csv(csv_filepath, index=False)
