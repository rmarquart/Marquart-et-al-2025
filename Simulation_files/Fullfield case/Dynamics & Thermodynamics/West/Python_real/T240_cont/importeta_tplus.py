# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 11:49:58 2023

@author: Ma-re Admin
"""

import os
import numpy as np

# Load the exported variable
num = str(os.environ.get("num"))

# Set the base path and 'numbers'
# basePath = os.getcwd()
basePath = "/scratch/mrqrut001/simulation/postdoc/Python_real"

# Create the folder path
folderPath = os.path.join(basePath, f'T{num}')
os.chdir(folderPath)

def is_numeric(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# Check if there are other numeric subfolders (including integers and floats)
subfolders = [folder for folder in os.listdir(folderPath) if is_numeric(folder) and folder != '0']
# print(subfolders)

if subfolders:
    highest_folder = max(subfolders, key=float)
    highest_folder_path = os.path.join(folderPath, highest_folder)

    # Check if 'eta' file exists in the highest folder
    h_file_path = os.path.join(highest_folder_path, 'eta')
    if os.path.exists(h_file_path):
        print(f"Opening 'eta' file in folder {highest_folder}:")
        with open(h_file_path, 'r') as fileID:
            for j in range(23):
                Headers = fileID.readline()

            # Moved this block outside of the loop to read data once
            Data = np.genfromtxt(fileID, dtype=float, max_rows=254016)
        
        # Moved this block outside of the loop to process and save the data once
        eta = Data
        numbers = folderPath[-3:]
        os.chdir('..')
        os.chdir('Storage eta/')
        os.chdir(f'T{num}')
        filename = f'eta{numbers}'
        np.save(filename, eta)
        os.chdir('..')
    else:
        print(f"'eta' file not found in folder {highest_folder}.")
else:
    print("No numeric subfolders found within the folder.")
