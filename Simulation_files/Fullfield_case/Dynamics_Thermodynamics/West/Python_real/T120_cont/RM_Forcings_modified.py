# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 12:53:25 2023

@author: Ma-re Admin
"""

#reading forcing data - high res (24h) ERA5 data only
import pandas as pd
import numpy as np

filename = 'ERA5_2022.xlsm'

# Read the Excel file
df1 = pd.read_excel(filename, sheet_name='ERA5')

# Extract the columns
Cl = df1['Cl'].values
qa = df1['qa'].values
qs = df1['qs'].values
Fsd_cloud = df1['Fsd_cloud'].values
U = df1['U'].values
V = df1['V'].values
P_rate = df1['P_rate'].values
Ta = df1['Ta'].values

# Read the Excel file
df2 = pd.read_excel(filename, sheet_name='OHF-BacktrackDrift')

# Read time series of ocean forcings
OHF = df2['OHF'].values
Sw = df2['sal'].values
Tfreez = df2['Tffr'].values
txtData = df2['Time']

fmt = '%Y-%m-%dT%H:%M:%S'

# Convert date-time strings to datetime objects
DateTime = [dt.timestamp() for dt in txtData]

# Enter your time step
# OPTIONS: 3600s, 2700s, 1800s, 1200s, 900s, 450s, 300s
t_orig = 3600
dt = 450
factor = dt / t_orig
k = dt / 3600

mult = 1 / k
divisions = round(len(DateTime) * mult - (mult - 1))

indices_to_remove = np.r_[0:33355, 52707:70073] #2022 dt=450

xi = np.linspace(min(DateTime), max(DateTime), divisions)
xi = np.delete(xi, indices_to_remove)

Cl = np.interp(xi, DateTime, Cl)
qa = np.interp(xi, DateTime, qa)
qs = np.interp(xi, DateTime, qs)
Fsd_cloud = np.interp(xi, DateTime, Fsd_cloud)
U = np.interp(xi, DateTime, U)
V = np.interp(xi, DateTime, V)
P_rate = np.interp(xi, DateTime, P_rate)
Ta = np.interp(xi, DateTime, Ta)
OHF = np.interp(xi, DateTime, OHF)
Sw = np.interp(xi, DateTime, Sw)
N = np.interp(xi, DateTime, DateTime)

Fsd_cloud = (Fsd_cloud / dt) * factor
Ua = np.sqrt(U**2 + V**2)
