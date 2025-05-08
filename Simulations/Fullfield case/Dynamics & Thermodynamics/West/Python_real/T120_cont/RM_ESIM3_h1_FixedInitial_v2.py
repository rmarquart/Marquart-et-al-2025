# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 14:31:34 2023

@author: Ma-re Admin
"""
import RM_Forcings_modified
import numpy as np
import pickle
import os
import pandas as pd
import shutil

num = str(os.environ.get("num"))

# Define the path to the 'storage.csv' file
csv_filepath = "/scratch/mrqrut001/simulation/postdoc/Python_real/Storage csv/T{}/storage.csv".format(num)

# Read the CSV file into a DataFrame
df = pd.read_csv(csv_filepath)

# Access the values from the DataFrame
iceFloeCategoryAverages = np.array(df['iceFloeCategoryAverages'])

# Change directory to 'Previous index values/h1/T{num}' folder
os.chdir(f'/scratch/mrqrut001/simulation/postdoc/Python_real/Previous index values/2022/dt=450s/h1/T{num}')

for i, category in enumerate(iceFloeCategoryAverages):

    if np.isnan(category):

        print(f'entry of h1 in category {i}: 0')
        print(f'outcome of h1 in category {i}: 0')
        
    else:

        filename = f'prev_index_values_{i}.pkl'
        
        # Print the name of the file that will be saved
        print(f"Loading file: {filename}")

        with open(filename, 'rb') as f:
            loaded_variables = pickle.load(f)

        print(f"'index' variable for {filename}: {loaded_variables['index']}")

        hi = iceFloeCategoryAverages[i]
        print(f'entry of h1 in category {i}: ', "{:.6f}".format(hi))
        
        #Station set-up
        nmax = 1                                 # maximum number of time steps (-)
        hs_prec_min = 0.02                       # minimum snow precipitation (m) to activate the model
        hs_min = 0.02                            # minimum snow thickness (m) to activate the model
        hmi_min = 0.05                           # minimum snow ice/superimposed ice thickness (m) to activate the model
        hi_min = 0.05                            # minimum sea ice thickness (m) to activate the model
        ZERO = 1e-5                              # ZERO OF THE MODEL
        Fwater = 11.7                            # oceanic heat fluxes (W m^-2)
        Tfreez = 271.28                          # seawater freezing temperature (K)

        #Parameters and constants
        alpha_ow = 0.15                          # seawater albedo (non-dim)
        alpha_melti = 0.55                       # intermediate layer albedo (non-dim)
        beta = 2.0                               # empir coeff from snow ice to ice (after Lepparanta, 1983) (non-dim)
        beta2 = 0.13                             # emp cost (after Untersteiner, 1961) (W m^-1)
        c0 = 2093.0                              # specific heat of fresh ice (J (kg*??deg)-1)
        c1 = 21.875                              # emp constant (K)
        c2 = 7.66                                # emp constant (K)
        c10 = 0.1                                # (m) (after Flato and Brown, 1996 for landfast sea ice of the Arctic)
        c11 = 0.44                               # (m^-0.28) (after Flato and Brown, 1996 for landfast sea ice of the Arctic)
        c12 = 0.075                              # (m^-2) (after Flato and Brown, 1996 for landfast sea ice of the Arctic)
        cail = 1.75e-3                           # bulk transfer coefficient (non-dim)
        cais = 1.20e-3                           # bulk transfer coefficient (non-dim)
        cpair = 1004.0                           # air heat capacity (J kg^-1 K^-1)
        cpw = 4186.0                             # specific heat of sea water (J kg^-1 K^-1)
        cs = 2093.0                              # specific heat of snow (J (kg K)-1)
        cwih = 2.8e-4                            # bulk heat trans coeff (after Omstedt and Wettlaufer, 1992) (non-dim)
        deltat = RM_Forcings_modified.dt         # time step (s)
        emis = 0.97                              # water/ice emissivity (non-dim)
        emp = 0.62197                            # empirical coefficient for specific humidity computation (mbar)
        gamm = 18.0                              # emp coeff (after Bitz and Lipscomb, 1999; otherwise 17.2 after Untersteiner, 1961) (J (K*kg)-1)
        infra_fr = 0.77                          # infrared partition as for water type II (after Jerlov, 1968) (non-dim)
        k0 = 2.03                                # fresh ice thermal conductivity (W m^-1 K^-1)
        ks_prec = 0.0641                         # precipit thermal conductivity (W m^-1 K^-1)
        kw = 0.563                               # sea water thermal conductivity (after Lepparanta, 1983) (W m^-1 K^-1)
        kmi = 0.9                                # snow ice thermal conductivity (after Lepparanta, 1983) (W m^-1 K^-1)
        ki = 1.90                                # Initial sea ice thermal conductivity (W m^-1 K^-1)
        k_ocean_red = 0.6667                     # seawater extinction coefficient (infrared, water type II, after Jerlov 1968) (m^-1)
        k_ocean_vis = 0.0714                     # seawater extinction coefficient (visible, water type II, after Jerlov 1968) (m^-1)
        h_mix = 10.0                             # mixed layer depth (m)
        Lai = 2.501e6                            # latent heat of vaporization (J kg^-1)
        L0i = 297000.0                           # latent heat of fusion of pure ice ()
        L0m = 315000.0                           # latent heat of fusion of meteoric ice ()
        L0s = 334000.0                           # latent heat of fusion of pure snow ()
        Lv = 2.501e6                             # latent heat of vaporization of fresh water at 273.15 K (J kg^-1)
        mu = 0.054                               # ratio between Tfr and brines salinity (after Assur, 1958) (??C)
        P0 = 1013.0                              # seawater pressure at the surface (mbar)
        q1 = 1.16378e7                           # emp coeff (kg m^-3)
        q2 = 5897.8                              # emp coeff (K)
        roa = 1.225                              # air density (kg m^-3)
        roi_bott = 0.9                           # pure sea ice density at the bottom (Kg m^-3)
        romi = 850.0                             # snow ice/superimposed ice density (kg m^-3)
        rosa = 1.68e-5                           # emp coeff (C??^-1*m*s^-1) (after Vancoppenolee et al, 2007; otherwise 1.68e-7 after Cox and Weeks, 1988)
        ros_prec = 250                           # snow precip density (kg m^-3)
        row = 1026.0                             # water density (kg m^-3)
        Sbr_end = 18.5185                        # pure brines density at the bottom (Kg m^-3)
        sigma = 5.68e-8                          # Stefan-Boltzmann constant (W m^-2 K)
        si_fract = 0.0                           # fraction of snow transformation in snow ice (non-dim) !!!IT CAN BE SET UP BY THE USER DEPENDING ON SNOW TYPE, ECC.)!!!
        ss_fract = 0.0                           # fraction of snow transformation in superimposed ice (after Cheng et al, 2006) (non-dim) !!!IT CAN BE SET UP BY THE USER DEPENDING ON SNOW TYPE, ECC.)!!!
        Ssnow = 0.0                              # snow salinity (per mill)
        Tb = 273.155                             # empirical temp (K)
        Tfr = Tfreez                             # freezing temp of sea ice (K)
        Tfrs = 273.15                            # freezing temp of snow/fresh water (K)
        Va = 0.015                               # gas volume fraction in the ice is fixed ()
        viola = 20.0                             # for salinity compuatation (after Vancoppenolee et al, 2007)
        vis_fr = 0.7                             # fraction of visible light penetrating the surface layer when SIM is on ()

        #Set the maximum number of iterations
        maxIterations = 100

        Tstore4 = np.zeros((3, maxIterations))
        Tstore2 = np.zeros((3, maxIterations))
              
        count_time_steps = 1                     # Keep track of time step

        Cl = RM_Forcings_modified.Cl[loaded_variables['index']]
        qa = RM_Forcings_modified.qa[loaded_variables['index']]
        qs = RM_Forcings_modified.qs[loaded_variables['index']]
        Fsd_cloud = RM_Forcings_modified.Fsd_cloud[loaded_variables['index']]
        U = RM_Forcings_modified.U[loaded_variables['index']]
        V = RM_Forcings_modified.V[loaded_variables['index']]
        P_rate = RM_Forcings_modified.P_rate[loaded_variables['index']]
        Ta = RM_Forcings_modified.Ta[loaded_variables['index']]
        OHF = RM_Forcings_modified.OHF[loaded_variables['index']]
        Sw = RM_Forcings_modified.Sw[loaded_variables['index']]

        Ua = np.sqrt(U**2 + V**2)

        #% MODEL START
        Sw = 34.22634  #Seawater salinity
                
        #Snow density
        if loaded_variables['T0'] < Tfr:
            ros = 300
            ros_new = 250
        else:
            ros = 350
            ros_new = 300

        #Snow precipitation
        if Ta < Tfr:
            if hi > hi_min:
                hs_prec = (P_rate / ros_prec) * deltat
            else:
                hs_prec = 0.0
        else:
            hs_prec = 0.0

        #Precipitation in the bucket model
        if loaded_variables['hs_prec_bucket'] > 0.01:
            delta_bucket = loaded_variables['hs_prec_bucket'] - ZERO
            hs_prec_bucket = hs_prec + ZERO
        elif loaded_variables['hs_prec_bucket'] <= 0.01:
            delta_bucket = 0.0
            hs_prec_bucket = loaded_variables['hs_prec_bucket'] + hs_prec

        ros_av = (loaded_variables['hs'] * ros + hs_prec_bucket * ros_new) / (loaded_variables['hs'] + hs_prec_bucket)

        #Thermal conductivity of new/old snow
        ks_new = (ros_new * 10**-3) ** 2 * 2.85  # (after Abel, 1892)
        ks = (ros * 10**-3) ** 2 * 2.85  # (after Abel, 1892)

        ks_av = (ks * loaded_variables['hs'] + ks_new * hs_prec_bucket) / (loaded_variables['hs'] + hs_prec_bucket)

        #Albedo computation (after Flato and Brown, 1996 for landfast sea ice of the Arctic)
        if loaded_variables['T0'] < Tfr:
            alpha_i = max(alpha_ow, c11 * hi**0.28 + 0.08)
            alpha_s = 0.75
            alpha_mi = 0.70  # (after Perovich, 1996 for compacted snow)
        else:
            alpha_i = min(alpha_melti, c12 * hi**2 + alpha_ow)
            alpha_s = 0.65
            alpha_mi = 0.56  # (after Perovich, 1996 for melting white ice)

        if (loaded_variables['hs'] + hs_prec_bucket) > 0.1:
            alpha = alpha_s
        elif loaded_variables['hmi'] > hmi_min:
            alpha = alpha_mi
        elif hi > hi_min:
            alpha = min(alpha_s, alpha_i + (loaded_variables['hs'] + hs_prec_bucket) * (alpha_s - alpha_i) / c10)
        else:
            alpha = alpha_ow
            
        #Extinction coefficient computation
        if loaded_variables['T0'] < Tfr:
            ks_daily_severeear = 25
            ks_cloudy = 25
            ks_snow = 25
                
            kmi_clear = 21.05
            kmi_cloudy = 17.75
            kmi_av = -3.3 * Cl + 21.05
            
            ksi_10_clear = 17.1
            ksi_10_cloudy = 10.5
            ksi_10_av = -6.6 * Cl + 17.1
            
            ksi_clear = 1.6
            ksi_cloudy = 1.5
            ksi_av = -0.1 * Cl + 1.6
            
        elif loaded_variables['T0'] >= Tfr:
            ks_clear = 15
            ks_cloudy = 15
            ks_snow = 15
            
            kmi_clear = 11.7
            kmi_cloudy = 9.8
            kmi_av = -1.9 * Cl + 11.7
            
            ksi_10_clear = 8.4
            ksi_10_cloudy = 4.6
            ksi_10_av = -3.8 * Cl + 8.4
            
            ksi_clear = 1.5
            ksi_cloudy = 1.4
            ksi_av = -0.1 * Cl + 1.5
            
        else:
            print('ERROR: EXTINCTION COEFFICIENT')

        #Heat fluxes not dependent on T0
        #Short wave
        Fs = (1 - alpha) * Fsd_cloud

        #Latent heat
        ea = qa * (P0 / emp)
        es = qs * (P0 / emp)

        Fla = cail * roa * Lai * Ua * (qs - qa)

        #Downward long wave
        if hi > hi_min:
            Fld = ((sigma * (Ta**4)) - 85.6) * (1 + 0.26 * Cl)  # (after Guest, 1997 for Antarctica)
        else:
            Fld = (sigma * Ta**4 * (0.653 - 0.00535 * ea)) * (1 + 0.1762 * Cl**2)  # (after Bignami et al., 1995 for Mediterranean Sea)
            Fld = (sigma * Ta**4 * (0.685 + 0.00452 * ea)) * (1 + 0.36 * Cl**2)  # (after Zapadka et al., 2007 for the Baltic)

        #Oceanic fluxes
        Fw = Fwater

        hs_tot = loaded_variables['hs']+hs_prec_bucket

        #%SURFACE TEMPERATURE ITERATION AND FLUXES     
        T0_star = np.zeros(maxIterations)
        Ts_star = np.zeros(maxIterations)
        Ti_star = np.zeros(maxIterations)
        Tsi_star = np.zeros(maxIterations)
            
        T = np.zeros(maxIterations)
            
        T0_star[0] = loaded_variables['T0']
        Ts_star[0] = loaded_variables['Ts']
        Tsi_star[0] = loaded_variables['Tsi']
        Ti_star[0] = loaded_variables['Ti']
            
        kice_surf = k0 + beta2 * loaded_variables['Si'] / (loaded_variables['Ti'] - 273.15)
        kice_bott = k0 + beta2 * loaded_variables['Sice_bott'] / (Tfr - 273.15)
        kice_bio = k0 + beta2 * loaded_variables['Sice_5'] / loaded_variables['Tice_5']

        ksnow_interm = (2 * (kmi * ks_av)) / (hs_tot * kmi + loaded_variables['hmi'] * ks_av)
        kinterm_ice = (2 * (ki * kmi)) / (loaded_variables['hmi'] * ki + hi * kmi)
        ksnow_ice = (2 * (ki * ks_av)) / (hs_tot * ki + hi * ks_av)

        Kice_bott = 2 * kice_bott / hi
        Kice_surf = 2 * kice_surf / hi
        Kice_bio = 2 * kice_bio / loaded_variables['hi_bio']
        Kbott_bio = 2 * kice_bott / loaded_variables['hi_bio']

        Kmi = 2 * kmi / loaded_variables['hmi']
        Ks = 2 * ks_av / hs_tot

        Ksnow_interm = ksnow_interm
        Kinterm_ice = kinterm_ice
        Ksnow_ice = ksnow_ice

        mus = deltat / (c0 * (ros_new * loaded_variables['hs_prec_bucket'] + ros * loaded_variables['hs']))
        mumi = deltat / (romi * c0 * loaded_variables['hmi'])
        cpi_bio = c0 + (L0s * mu * loaded_variables['Sice_bio']) / (loaded_variables['Tice_bio']**2)
        mui_bio = deltat / (loaded_variables['ro_sice_bulk'] * 10**3 * cpi_bio * loaded_variables['hi_bio'])

        #% ITERATION START
            
        #Set the convergence tolerance
        tolerance = 0.01
            
        for j in range(1, maxIterations):
            Flu_it = emis * sigma * T0_star[j-1]**4
            Fse_it = roa * cpair * cais * Ua * (T0_star[j-1] - Ta)
            Fl_it = Flu_it - Fld
            F_it = - Fs + Fl_it + Fla + Fse_it
            
            la = 4 * emis * sigma * T0_star[j-1]**3
            sen = roa * cpair * cais * Ua
            lat = cail * roa * Lai * Ua * ((c1 * (Tb - c2) * qs) / ((T0_star[j-1] + Tb - c2)**2))

            cons = ks_av / (hs_tot / 2)
            coni = kice_surf / (hi / 2)
            conmi = kmi / (loaded_variables['hmi'] / 2)
                
            #if i == 1 and j == 1:
            #    Ti_sqr = (Ti_star[j-1] - 273.15) * (Ti_star[j-1] - 273.15)
            #else:
            #    Ti_sqr = (Ti_star_array[np.nonzero(Ti_star_array)][-1] - 273.15) * (Ti_star_array[np.nonzero(Ti_star_array)][-2] - 273.15)
            
            #cpi_it = c0 + ((L0s * mu * Si[i - 1]) / (Ti_sqr))
            cpi_it = c0 + ((L0s * mu * loaded_variables['Si']) / ((Ti_star[j-1] - 273.15) * (Ti_star[0] - 273.15)))
                
            mui = deltat / (loaded_variables['ro_sice_bulk'] * 10**3 * cpi_it * hi)

            if (loaded_variables['hs'] + hs_prec_bucket) > hs_min:
                
                #SEA ICE MODEL ON/SLAB OCEAN OFF
                if loaded_variables['hmi'] > hmi_min:
                    strI = '1'
                
                elif loaded_variables['hmi'] <= hmi_min:
                    strI = '2'
                    
                    I0_it = vis_fr * Fs * np.exp(-ks_snow * (loaded_variables['hs'] + hs_prec_bucket))
                    IM_it = I0_it
                    
                    if hi >= 0.1:
                        ISI_10_it = IM_it * np.exp(-ksi_10_av * 0.1)
                        ISI_it = ISI_10_it * np.exp(-ksi_av * (hi - 0.1))
                    elif hi < 0.1:
                        ISI_10_it = IM_it * np.exp(-ksi_10_av * hi)
                        ISI_it = ISI_10_it

                    A = np.array([
                        [(1 - Ks / (la + sen + lat + cons)), (Ks / (la + sen + lat + cons)), 0],
                        [(mus * Ks), (1 - mus * (Ksnow_ice + Ks)), (mus * Ksnow_ice)],
                        [0, (mui * Ksnow_ice), (1 - mui * (Ksnow_ice + Kice_bott))]
                    ])

                    b = np.array([
                        ((-F_it - I0_it) / (la + sen + lat + cons)),
                        (I0_it * mus),
                        (mui * Kice_bott * Tfr + mui * ISI_it)])

                    Tprev = np.array([T0_star[j-1], Ts_star[0], Ti_star[0]])

                    T = np.dot(A, Tprev) + b
                                    
                    Tstore2[:,j] = T 
                                                                                                                  
                    T0_star[j] = Tstore2[0,j]
                    Ts_star[j] = Tstore2[1,j]
                    Ti_star[j] = Tstore2[2,j]

                    T0_iter = T0_star[np.nonzero(T0_star)][-1]
                    Ts_iter = Ts_star[np.nonzero(Ts_star)][-1]
                    Ti_iter = Ti_star[np.nonzero(Ti_star)][-1]
                                                            
                    Tmix_star = Tfr

                    if Ts_iter > Tfr:
                        Ts_iter = Tfr

                    if Ti_iter > Tfr:
                        Ti_iter = Tfr
                        Ts_iter = Tfr
                        
                    if T0_iter > Tfr:
                        T0_iter = Tfr
                            
                        B = np.array([
                            [0, 0, 0],
                            [(mus * Ks), (1 - mus * (Ksnow_ice + Ks)), (mus * Ksnow_ice)],
                            [0, (mui * Ksnow_ice), (1 - mui * (Ksnow_ice + Kice_bott))]
                        ])

                        c = np.array([0, (I0_it * mus), (mui * Kice_bott * Tfr + mui * ISI_it)])
                            
                        Tprev_melt = np.array([Tfr, Ts_star[0], Ti_star[0]])  
                            
                        Tmelt = np.dot(B, Tprev_melt) + c
                                                                   
                        Ts_iter = Tmelt[1]
                            
                        if Ts_iter < Tfr:
                            Ti_iter = Tmelt[2]
                        else:
                            Ts_iter = Tfr
                            Ti_iter = (mui * Ksnow_ice) * Tfr + (1 - mui * (Ksnow_ice + Kice_bott)) * loaded_variables['Ti'] + (mui * Kice_bott * Tfr + mui * ISI_it)
                            if Ti_iter > Tfr:
                                Ti_iter = Tfr
                        
                    Tsi_iter=Tfr

            elif (loaded_variables['hs'] + hs_prec_bucket) <= hs_min:
                if loaded_variables['hmi'] > hmi_min:
                    strI = '3'
                
                elif loaded_variables['hmi'] <= hmi_min:
                    if hi > hi_min:
                        strI = '4'
                            
                        if hi >= 0.1:
                            ISI_10_it = vis_fr * Fs * np.exp(-ksi_10_av * 0.1)
                            ISI_it = ISI_10_it * np.exp(-ksi_av * (hi - 0.1))
                        elif hi < 0.1:
                            ISI_10_it = vis_fr * Fs * np.exp(-ksi_10_av * hi)
                            ISI_it = ISI_10_it
                                 
                        I0_it = ISI_10_it
                        IM_it = ISI_10_it
                                 
                        A = np.array([
                            [(1 - Kice_surf / (la + sen + lat + coni)), (Kice_surf / (la + sen + lat + coni)), 0],
                            [(mui * Kice_surf), (1 - mui * (Kice_surf + Kice_bott)), 0],
                            [0, 0, 0]
                        ])
                                                     
                        b = np.array([((-F_it - I0_it) / (la + sen + lat + coni)), (mui * Kice_bott * Tfr + mui * ISI_it), 0])
                                 
                        Tprev = np.array([T0_star[j-1], Ti_star[0], Ti_star[0]])
                            
                        T = np.dot(A, Tprev) + b
                            
                        #if i in [8145,8146]:
                        #    rounded_values_2 = [A,b,Tprev]
                        #    print(rounded_values_2)
                                                
                        #if i == 8147:
                        #    breakpoint()
                            
                        Tstore4[:,j] = T
                                                     
                        T0_star[j] = Tstore4[0,j]
                        Ti_star[j] = Tstore4[1,j]                   
                                                                       
                        T0_iter = T0_star[np.nonzero(T0_star)][-1]
                        Ti_iter = Ti_star[np.nonzero(Ti_star)][-1]
                                 
                        Tmix_star = Tfr
                                 
                        if Ti_iter > Tfr:
                            Ti_iter = Tfr
                                     
                        if T0_iter > Tfr:
                            T0_iter = Tfr
                            Ti_iter = (mui * Kice_surf) * Tfr + (1 - mui * (Kice_surf + Kice_bott)) * loaded_variables['Ti'] + (mui * Kice_bott * Tfr + mui * ISI_it)
                            if Ti_iter > Tfr:
                                Ti_iter = Tfr
                                     
                        Ts_star[j] = Tfr
                        Ts_iter = Ts_star[np.nonzero(Ts_star)][-1]
                        Tsi_iter = Tfr

                    elif hi <= hi_min:
                        strI = '5'

                        I0_it = Fs * ((infra_fr * np.exp(-k_ocean_red * h_mix)) + ((1 - infra_fr) * np.exp(-k_ocean_vis * h_mix)))
                        Tmix_star = loaded_variables['Tmix'] - ((F_it + I0_it + loaded_variables['R']) / (h_mix * row * cpw)) * deltat
                        
                        T0_star[j] = Tmix_star
                            
                        Ts_star[j] = Tfr
                            
                        Tsi_star[j] = Tfr
                        Ti_star[j] = Tfr
                                                
                        T0_iter = T0_star[np.nonzero(T0_star)][-1]
                        Ts_iter = Ts_star[np.nonzero(Ts_star)][-1]
                        Tsi_iter = Tsi_star[np.nonzero(Tsi_star)][-1]
                        Ti_iter = Ti_star[np.nonzero(Ti_star)][-1]
                                               
                else:
                    print('ERROR ITERATION')
                    
            if abs(T0_star[j] - T0_star[j-1]) < tolerance:
                break #Go to the next time step
            
        if j > 20:
            print('WARNING: ITERATION > 20')

        if j == maxIterations:
            print('NOT CONVERGENT')
                                
        #if j == maxIterations:
            #break #Stop the sea ice model if not convergent
               
        T0 = T0_iter       # NEW SURFACE TEMPERATURE
        Ts = Ts_iter       # NEW SNOW TEMPERATURE
        Tsi = Tsi_iter     # NEW SNOW ICE/SUPERIMPOSED ICE TEMPERATURE
        Ti = Ti_iter       # NEW SEA ICE TEMPERATURE
        Tmix = Tmix_star   # NEW MIX LAYER TEMPERATURE
          
        if loaded_variables['hmi'] > hmi_min:
            Tice = ((kice_surf * Ti / (hi / 2) + ksnow_interm * Tsi / (loaded_variables['hmi'] / 2)) / (kice_surf / (hi / 2) + ksnow_interm / (loaded_variables['hmi'] / 2))) - 273.15
        elif loaded_variables['hmi'] <= hmi_min:
            if (loaded_variables['hs'] + hs_prec_bucket) > hs_min:
                Tice = ((kice_surf * Ti / (hi / 2) + ks_av * Ts / (hs_tot / 2)) / (kice_surf / (hi / 2) + ks_av / (hs_tot / 2))) - 273.15
            elif (loaded_variables['hs'] + hs_prec_bucket) <= hs_min:
                if hi > hi_min:
                    Tice = T0 - 273.1499
                elif hi <= hi_min:
                    Tice = Tfr - 273.15
            
        #Keep track of Ti_iter T0_iter Ts_iter            
        #empty_Ti_star[:,i] = np.transpose(Ti_star)
        #empty_T0_star[:,i] = np.transpose(T0_star)
        #empty_Ts_star[:,i] = np.transpose(Ts_star)

        #empty_Ti_star_trans = np.transpose(empty_Ti_star)
        #empty_T0_star_trans = np.transpose(empty_T0_star) 
        #empty_Ts_star_trans = np.transpose(empty_Ts_star)
            
        #Ti_star_array = empty_Ti_star_trans[i,:]

        #Check if variables exist     
        variables_to_delete = ['T0_star', 'Ts_star', 'Ti_star', 'Tsi_star', 'Tmix_star', 'Tprev', 'T', 'A', 'a', 'B', 'b', 'C', 'c', 'D', 'd', 'e', 'Flu_it', 'Fse_it', 'Fl_it', 'F_it', 'la', 'sen', 'lat', 'ISI_it', 'ISI_10_it', 'IMI_it', 'I0_it', 'Tprev_melt', 'Tmelt', 'Tprev_melt_2', 'Tmelt_2', 'T0_iter', 'Ts_iter', 'Tsi_iter', 'Ti_iter', 'cons', 'coni', 'conmi', 'mui', 'cpi_it']

        for variable_name in variables_to_delete:
            if variable_name in locals():
                del locals()[variable_name]
                    
        #Recompute surface fluxes
        Flu = emis * sigma * T0**4  # UPWARD LONGWAVE RADIATION
        Fl = Flu - Fld  # NET LONGWAVE RADIATION
        Fse = roa * cpair * cais * Ua * (T0 - Ta)  # SENSIBLE HE
        F = -Fs + Fl + Fla + Fse  # NET SURFACE FLUXES
        print(f'Net surface flux: {F:.6f}')
        R = 0.0  # INITIAL EXTRA HEAT IS SET TO ZERO

        Qs = ros * ((c0 * (Tfrs - Ts)) + L0s)
        Qmi = romi * ((c0 * (Tfrs - Tsi)) + L0s)
        Qi_surf = loaded_variables['ro_sice_surf'] * 10**3 * ((c0 * (Tfr - Ti)) + L0s * (1 + mu * loaded_variables['Si'] / (Ti - 273.15)))
        Qi_bott = loaded_variables['ro_sice_bott'] * 10**3 * ((c0 * (Tfr - Ti)) + L0s * (1 + mu * loaded_variables['Si'] / (Ti - 273.15)))

        #%Sea ice growth/decay
        if T0 < Tfr:
            if loaded_variables['hmi'] > hmi_min:
                
                if (loaded_variables['hs'] + hs_prec_bucket) > hs_min:
                    strII = '1'
                elif (loaded_variables['hs'] + hs_prec_bucket) <= hs_min:
                    strII = '2'
                
            elif loaded_variables['hmi'] <= hmi_min:
                    
                if (loaded_variables['hs'] + hs_prec_bucket) > hs_min:
                    strII = '3'
                    Fsurf = Ks * (Ts - T0)
                    Fbott = Kice_bott * (Tfr - Ti)
                        
                    I0 = vis_fr * Fs * np.exp(-ks_snow * (loaded_variables['hs'] + hs_prec_bucket))
                    IM = I0
                        
                    if hi >= 0.1:
                        ISI_10 = IM * np.exp(-ksi_10_av * 0.1)
                        ISI = ISI_10 * np.exp(-ksi_av * (hi - 0.1))
                    elif hi < 0.1:
                        ISI_10 = IM * np.exp(-ksi_10_av * hi)
                        ISI = ISI_10
                        
                    delta_sublim = (Fla * deltat) / (ros_new * Lv - Qs)
                    deltahs_melt_surf = 0.0
                    deltahi_melt_surf = 0.0
                    deltahmi_melt_surf = 0.0
                    deltahi_bott = deltat * (Fbott - Fw - ISI) / Qi_bott
                    
                    if (loaded_variables['hs'] + hs_prec_bucket - delta_sublim) > (hi * (row - loaded_variables['ro_sice_surf'] * 10**3) / ros_av - hs_prec_bucket * ros_new / ros):
                        hs = loaded_variables['hs'] + delta_bucket - delta_sublim - beta * ((-row * hi + ros * loaded_variables['hs'] + 10**3 * loaded_variables['ro_sice_surf'] * hi + ros_new * hs_prec_bucket) / (-romi + row + beta * ros_av))
                        hmi_new = ((-row * hi + ros * (loaded_variables['hs'] - delta_sublim) + 10**3 * loaded_variables['ro_sice_surf'] * hi + ros_new * hs_prec_bucket) / (-romi + row + beta * ros_av)) * si_fract
                        hmi = loaded_variables['hmi'] + hmi_new
                    else:
                        hs = loaded_variables['hs'] + delta_bucket - delta_sublim
                        hmi_new = 0.0
                        hmi = loaded_variables['hmi']
                            
                    if hs < ZERO:
                        hs = hs_prec_bucket
                        hs_prec_bucket = 0.0
                        
                    hi = hi + deltahi_bott
                        
                elif (loaded_variables['hs'] + hs_prec_bucket) <= hs_min:
                    strII = '4'
                    Fsurf = Kice_surf * (Ti - T0)
                    Fbott = Kice_bott * (Tfr - Ti)
                        
                    if hi >= 0.1:
                        ISI_10 = vis_fr * Fs * np.exp(-ksi_10_av * 0.1)
                        ISI = ISI_10 * np.exp(-ksi_av * (hi - 0.1))
                    elif hi < 0.1:
                        ISI_10 = vis_fr * Fs * np.exp(-ksi_10_av * hi)
                        ISI = ISI_10
                        
                    I0 = ISI_10
                    IM = ISI_10
                        
                    if hi > hi_min:
                        delta_sublim = (Fla * deltat) / (loaded_variables['ro_sice_surf'] * Lv - Qi_surf)
                        deltahi_bott = deltat * (Fbott - Fw - ISI) / Qi_bott
                    else:
                        delta_sublim = 0.0
                        ISI = Fs * ((infra_fr * np.exp(-k_ocean_red * h_mix)) + ((1 - infra_fr) * np.exp(-k_ocean_vis * h_mix)))
                        deltahi_bott = deltat * (F + ISI) / Qi_bott
                        
                    deltahs_melt_surf = 0.0
                    deltahmi_melt_surf = 0.0
                    deltahi_melt_surf = 0.0
                    
                    hs = loaded_variables['hs'] + delta_bucket
                    
                    hmi_new = 0.0
                    hmi = loaded_variables['hmi']
                    hi = hi + deltahi_bott - delta_sublim
                    
                    if hi < ZERO:
                        R = (loaded_variables['ro_sice_bott'] * 10**3 * ((c0 * (Tfr - Ti)) + L0s * (1 + loaded_variables['Si'] / (Ti - 273.15)))) * (hi - ZERO) / deltat
                        hi = ZERO
                        hs = ZERO
                        hmi = ZERO
                        
            else:
                print('ERROR GROWTH SEASON')

            #% HALODYNAMIC SUBMODEL FOR THE GROWTH SEASON
            if hi > ZERO:
                vi = deltahi_bott / deltat
                keff = 0.12 / (0.12 + 0.88 * np.exp(-4.2 * 10**4 * vi))  # with 4.2*10^4 in s/cm for normal salinity sea water >30, it is usually used the parameterization by Nawako and Shina, 1983
                Sice_bott = keff * Sw
                Sbr_bott = -(Tfr - 273.1499) / mu
                Vbr_bott = Sice_bott / Sbr_bott
                Q = 0.0
                    
                if hmi_new > 0.0:
                    snow_fr = beta / (loaded_variables['ro_sice_surf'] * 10**3 / ros_av)
                    sea_water_fr = 1 - snow_fr
                    Ssnowice_new = (hmi_new * sea_water_fr * Sw) / hmi_new
                    Ssnowice = (loaded_variables['Ssnowice'] * loaded_variables['hmi'] + Ssnowice_new * hmi_new) / hmi
                else:
                    snow_fr = loaded_variables['snow_fr']
                    sea_water_fr = loaded_variables['sea_water_fr']
                    Ssnowice_new = 0.0
                    Ssnowice = loaded_variables['Ssnowice']
                    
                if loaded_variables['Vbr_ice'] > 0.05:
                    Sice_mix = ((hi - hi) * Sice_bott + hi * loaded_variables['Sice']) / hi
                    Sice = Sice_mix + (rosa * (1 - viola * loaded_variables['Vbr_ice']) * (Tfr - 273.15 - Tice)) * deltat
                    Sbr_ice = -Tice / mu
                    Vbr_ice = Sice / Sbr_ice
                    Sice_5 = Sice
                    Sbr_5 = Sbr_ice
                    Tice_5 = Tice
                    Tice_bio = (loaded_variables['ki_bio_bott'] * (Tfr - 273.1499) + loaded_variables['ki_ice_bio'] * Tice) / (loaded_variables['ki_bio_bott'] + loaded_variables['ki_ice_bott'])
                
                elif loaded_variables['Vbr_ice'] <= 0.05:
                    Sice_mix = ((hi - hi) * Sice_bott + hi * loaded_variables['Sice']) / hi
                    Sice = Sice_mix
                    Sbr_ice = -Tice / mu
                    Vbr_ice = Sice / Sbr_ice
                    Sice_5_mix = (loaded_variables['hi_5'] * loaded_variables['Sice_5'] + (hi - hi) * Sice_bott) / (hi - hi + loaded_variables['hi_5'])
                    Sice_5 = Sice_5_mix + (rosa * (1 - viola * 0.0499) * (Tfr - 273.15 - loaded_variables['Tice_5'])) * deltat
                    Sbr_5 = Sice_5 / 0.0499
                    Tice_5 = -mu * Sbr_5
                    Tice_bio = (loaded_variables['ki_bio_bott'] * (Tfr - 273.1499) + loaded_variables['ki_ice_5'] * Tice_5) / (loaded_variables['ki_bio_bott'] + loaded_variables['ki_ice_5'])
                
                hi_5 = (hi * loaded_variables['ki_5_bott'] * (Tice_5 - (Tfr - 273.1499))) / ((loaded_variables['ki_ice_5'] * (Tice - Tice_5)) + (loaded_variables['ki_5_bott'] * (Tice_5 - (Tfr - 273.1499))))
                hi_bio = hi_5
                Sice_bio_mix = ((hi - hi) * Sice_bott + loaded_variables['hi_bio'] * loaded_variables['Sice_bio']) / (hi - hi + loaded_variables['hi_bio'])
                Sice_bio = Sice_bio_mix + (rosa * (1 - viola * loaded_variables['Vbr_bio']) * (Tfr - 273.15 - Tice_bio)) * deltat
                Sbr_bio = -Tice_bio / mu
                Vbr_bio = Sice_bio / Sbr_bio
                
                Si_mix = ((hi - hi) * Sice_bott + (hi / 2) * loaded_variables['Si']) / (hi - hi + hi / 2)

                if loaded_variables['Vbr_i'] <= 0.05:
                    Si = Si_mix
                else:
                    Si = Si_mix + (rosa * (1 - viola * loaded_variables['Vbr_i']) * (273.1499 - Ti)) * deltat
                    
                Sbr_i = -(Ti - 273.1499) / mu
                Vbr_i = Si / Sbr_i

            elif hi <= ZERO:
                snow_fr = 0.0
                sea_water_fr = 0.0
                Ssnowice_new = 0.0
                Ssnowice = 0.0
                
                Q = ZERO
                Sice = Sw
                Sice_bott = Sw
                
                Sice_mix = Sw
                Sice_5_mix = Sw
                Sice_bio_mix = Sw
                
                Sbr_ice = Sw
                Vbr_ice = 1.0
                
                Sbr_5 = Sw
                Sice_5 = Sw
                Tice_5 = Tfr - 273.15
                hi_5 = ZERO
                    
                Si = Sw
                Vbr_i = 1.0
                Sbr_i = Sw
                    
                Sbr_bio = Sw
                Sice_bio = Sw
                Tice_bio = Tfr - 273.15
                Vbr_bio = 1.0
                hi_bio = hi_5
                ISI_bio = 0.0
                
                Sice_bott = Sw
                Sbr_bott = Sw
                Vbr_bott = 1.0

        elif T0 >= Tfr:
            keff = 1.000
            vi = 0.0
            hmi_new = 0.0
            snow_fr = 0.0
            sea_water_fr = 0.0
            Ssnowice_new = 0.0
            
            if (loaded_variables['hs'] + hs_prec_bucket) > hs_min:
                strIII = '1'
                T0 = Tfr
                I0 = vis_fr * Fs * np.exp(-ks_snow * (loaded_variables['hs'] + hs_prec_bucket))
                Fsurf = Ks * (Ts - T0)
                Fbott = Kice_bott * (Tfr - Ti)
                
                if loaded_variables['hmi'] > hmi_min:
                    IM = I0 * np.exp(-kmi_av * loaded_variables['hmi'])
                elif loaded_variables['hmi'] <= hmi_min:
                    hmi = ZERO
                    IM = I0
                
                if hi >= 0.1:
                    ISI_10 = IM * np.exp(-ksi_10_av * 0.1)
                    ISI = ISI_10 * np.exp(-ksi_av * (hi - 0.1))
                elif hi < 0.1:
                    ISI_10 = IM * np.exp(-ksi_10_av * hi)
                    ISI = ISI_10
                    
                if Fsurf >= (F + I0):
                    deltahs_melt_surf = deltat * (F + I0 - Fsurf) / Qs
                else:
                    deltahs_melt_surf = 0.0
                    
                delta_sublim = (Fla * deltat) / (ros_av * Lv - Qs)
                deltahi_melt_surf = 0.0
                deltahmi_melt_surf = 0.0
                deltahi_bott = deltat * (Fbott - Fw - ISI) / Qi_bott
                
                hs = loaded_variables['hs'] + deltahs_melt_surf + delta_bucket - delta_sublim
                hmi = loaded_variables['hmi'] - deltahs_melt_surf * ros / romi * ss_fract
                hi = hi + deltahi_bott
                
                if hs < hs_min:
                    hs = hs_prec_bucket + (hs - ZERO)
                    hs_prec_bucket = 0.0
                    
                if hmi > hmi_min:
                    Ssnowice = (-deltahs_melt_surf * 0.20 * Ssnow + hmi * loaded_variables['Ssnowice']) / (hmi - deltahs_melt_surf)
                else:
                    Ssnowice = 0.0
                    
                if loaded_variables['Vbr_ice'] >= 0.05:
                    Q = 0.20 * ros * (-deltahs_melt_surf) / deltat
                    Sbr_ice_star = loaded_variables['Sbr_ice'] - ((Q / (loaded_variables['Vbr_ice'] * loaded_variables['ro_br'] * 10**3)) * (loaded_variables['Sbr_ice'] - Sbr_end)) * deltat
                    Sice = loaded_variables['Sice'] + (Sbr_ice_star - loaded_variables['Sbr_ice']) * loaded_variables['Vbr_ice']
                    Vbr_ice = -mu * Sice / Tice
                    Sbr_ice = Sice / Vbr_ice
                    Sbr_i_star = loaded_variables['Sbr_i'] - ((Q / (loaded_variables['Vbr_i'] * loaded_variables['ro_br_i'] * 10**3)) * (loaded_variables['Sbr_i'] - Sbr_end)) * deltat
                    Si = loaded_variables['Si'] + (Sbr_i_star - loaded_variables['Sbr_i']) * loaded_variables['Vbr_i']
                    Vbr_i = -mu * Si / (Ti - 273.15)
                    Sbr_i = Si / Vbr_i
                    Sbr_5_star = loaded_variables['Sbr_5'] - ((Q / (0.05 * loaded_variables['ro_br_5'] * 10**3)) * (loaded_variables['Sbr_5'] - Sbr_end)) * deltat
                    Sice_5 = loaded_variables['Sice_5'] + (Sbr_5_star - loaded_variables['Sbr_5']) * 0.05
                    Tice_5 = -mu * Sice_5 / 0.05
                    Sbr_5 = Sice_5 / 0.05
                    hi_5 = (hi * loaded_variables['ki_5_bott'] * (Tice_5 - (Tfr - 273.15))) / ((loaded_variables['ki_ice_5'] * (Tice - Tice_5)) + (loaded_variables['ki_5_bott'] * (Tice_5 - (Tfr - 273.15))))
                    Sbr_bio_star = loaded_variables['Sbr_bio'] - ((Q / (loaded_variables['Vbr_bio'] * loaded_variables['ro_br_bio'] * 10**3)) * (loaded_variables['Sbr_bio'] - Sbr_end)) * deltat  # DESALINATION IN THE ACTIVE BIOLOGICAL SYSTEM
                    Sice_bio = loaded_variables['Sice_bio'] + (Sbr_bio_star - loaded_variables['Sbr_bio']) * loaded_variables['Vbr_bio']
                    Tice_bio = (loaded_variables['ki_bio_bott'] * (Tfr - 273.15) + loaded_variables['ki_ice_bio'] * Tice) / (loaded_variables['ki_bio_bott'] + loaded_variables['ki_ice_bott'])
                    Vbr_bio = -mu * Sice_bio / Tice_bio
                    Sbr_bio = Sice_bio / Vbr_bio
                    hi_bio = hi_5
                    Sbr_bott_star = loaded_variables['Sbr_bott'] - ((Q / (loaded_variables['Vbr_bott'] * loaded_variables['ro_br_bott'] * 10**3)) * (loaded_variables['Sbr_bott'] - Sbr_end)) * deltat  # DESALINATION IN THE BOTTOM LAYER
                    Sice_bott = loaded_variables['Sice_bott'] + (Sbr_bott_star - loaded_variables['Sbr_bott']) * loaded_variables['Vbr_bott']
                    Vbr_bott = -mu * Sice_bott / (Tfr - 273.1499)
                    Sbr_bott = Sice_bott / Vbr_bott
                    
                elif loaded_variables['Vbr_ice'] < 0.05:
                    Q = 0.0
                    Sice = loaded_variables['Sice']
                    Sbr_ice = -Tice / mu
                    Vbr_ice = Sice / Sbr_ice
                    Si = loaded_variables['Si']
                    Sbr_i = -(Ti - 273.15) / mu
                    Vbr_i = Si / Sbr_i
                    Sice_5 = loaded_variables['Sice_5']
                    Sbr_5 = Sice_5 / 0.05
                    Tice_5 = -mu * Sbr_5
                    hi_5 = (hi * loaded_variables['ki_5_bott'] * (Tice_5 - (Tfr - 273.15))) / ((loaded_variables['ki_ice_5'] * (Tice - Tice_5)) + (loaded_variables['ki_5_bott'] * (Tice_5 - (Tfr - 273.15))))
                    Tice_bio = (loaded_variables['ki_bio_bott'] * (Tfr - 273.15) + loaded_variables['ki_ice_5'] * Tice_5) / (loaded_variables['ki_bio_bott'] + loaded_variables['ki_ice_5'])
                    Sice_bio = loaded_variables['Sice_bio']
                    Sbr_bio = -Tice_bio / mu
                    Vbr_bio = Sice_bio / Sbr_bio
                    hi_bio = hi_5
                    Sice_bott = loaded_variables['Sice_bott']
                    Sbr_bott = -(Tfr - 273.1499) / mu
                    Vbr_bott = Sice_bott / Sbr_bott
                    
            elif (loaded_variables['hs'] + hs_prec_bucket) <= hs_min:
                
                if loaded_variables['hmi'] > hmi_min:
                    strIII = '2'
                        
                elif loaded_variables['hmi'] <= hmi_min:
                    if hi > ZERO:
                        strIII = '3'
                        hs = ZERO
                        hmi = ZERO
                        T0 = Tfr
                            
                        if hi >= 0.1:
                            ISI_10 = vis_fr * Fs * np.exp(-ksi_10_av * 0.1)
                            ISI = ISI_10 * np.exp(-ksi_av * (hi - 0.1))
                        elif hi < 0.1:
                            ISI_10 = vis_fr * Fs * np.exp(-ksi_10_av * hi)
                            ISI = ISI_10
                        
                        I0 = ISI_10
                        IM = I0
                        
                        Fsurf = Kice_surf * (Ti - T0)
                        Fbott = Kice_bott * (Tfr - Ti)
                        
                        delta_sublim = (Fla * deltat) / (loaded_variables['ro_sice_surf'] * Lv - Qi_surf)
                        deltahi_bott = deltat * (Fbott - Fw - ISI) / Qi_bott
                    
                        if Fsurf >= (F + I0):
                            deltahi_melt_surf = deltat * (F + I0 - Fsurf) / Qi_surf
                        else:
                            deltahi_melt_surf = 0.0
                            
                        deltahs_melt_surf = 0.0
                        deltahmi_melt_surf = 0.0
                        
                        hs = ZERO + delta_bucket
                        hmi = ZERO
                        hi = hi + deltahi_melt_surf + deltahi_bott - delta_sublim
                        
                        if hi < ZERO:
                            R = (Qi_surf * (hi - ZERO)) / deltat
                            hi = ZERO
                            
                        Ssnowice = 0.0
                        
                        if hi > ZERO:
                            if loaded_variables['Vbr_ice'] >= 0.05:
                                Q = 0.20 * loaded_variables['ro_sice_surf'] * (-deltahi_melt_surf) / deltat
                                Sbr_ice_star = loaded_variables['Sbr_ice'] - ((Q / (loaded_variables['Vbr_ice'] * loaded_variables['ro_br'] * 10**3)) * (loaded_variables['Sbr_ice'] - Sbr_end)) * deltat
                                Sice = loaded_variables['Sice'] + (Sbr_ice_star - loaded_variables['Sbr_ice']) * loaded_variables['Vbr_ice']
                                Vbr_ice = -mu * Sice / Tice
                                Sbr_ice = Sice / Vbr_ice
                                Sbr_i_star = loaded_variables['Sbr_i'] - ((Q / (loaded_variables['Vbr_i'] * loaded_variables['ro_br_i'] * 10**3)) * (loaded_variables['Sbr_i'] - Sbr_end)) * deltat
                                Si = loaded_variables['Si'] + (Sbr_i_star - loaded_variables['Sbr_i']) * loaded_variables['Vbr_i']
                                Vbr_i = -mu * Si / (Ti - 273.15)
                                Sbr_i = Si / Vbr_i
                                Sbr_5_star = loaded_variables['Sbr_5'] - ((Q / (0.05 * loaded_variables['ro_br_5'] * 10**3)) * (loaded_variables['Sbr_5'] - Sbr_end)) * deltat
                                Sice_5 = loaded_variables['Sice_5'] + (Sbr_5_star - loaded_variables['Sbr_5']) * 0.05
                                Tice_5 = (loaded_variables['ki_5_bott'] * (Tfr - 273.15) + loaded_variables['ki_ice_5'] * Tice) / (loaded_variables['ki_5_bott'] + loaded_variables['ki_ice_5'])
                                Vbr_5 = -mu * Sice_5 / Tice_5
                                Sbr_5 = Sice_5 / Vbr_5
                                hi_5 = (hi * loaded_variables['ki_5_bott'] * (Tice_5 - (Tfr - 273.15))) / ((loaded_variables['ki_ice_5'] * (Tice - Tice_5)) + (loaded_variables['ki_5_bott'] * (Tice_5 - (Tfr - 273.15))))
                                Sbr_bio_star = loaded_variables['Sbr_bio'] - ((Q / (loaded_variables['Vbr_bio'] * loaded_variables['ro_br_bio'] * 10**3)) * (loaded_variables['Sbr_bio'] - Sbr_end)) * deltat
                                Sice_bio = loaded_variables['Sice_bio'] + (Sbr_bio_star - loaded_variables['Sbr_bio']) * loaded_variables['Vbr_bio']
                                Tice_bio = (loaded_variables['ki_bio_bott'] * (Tfr - 273.15) + loaded_variables['ki_ice_bio'] * Tice) / (loaded_variables['ki_bio_bott'] + loaded_variables['ki_ice_bott'])
                                Vbr_bio = -mu * Sice_bio / Tice_bio
                                Sbr_bio = Sice_bio / Vbr_bio
                                hi_bio = hi_5
                                Sbr_bott_star = loaded_variables['Sbr_bott'] - ((Q / (loaded_variables['Vbr_bott'] * loaded_variables['ro_br_bott'] * 10**3)) * (loaded_variables['Sbr_bott'] - Sbr_end)) * deltat
                                Sice_bott = loaded_variables['Sice_bott'] + (Sbr_bott_star - loaded_variables['Sbr_bott']) * loaded_variables['Vbr_bott']
                                Vbr_bott = -mu * Sice_bott / (Tfr - 273.1499)
                                Sbr_bott = Sice_bott / Vbr_bott
                            
                            elif loaded_variables['Vbr_ice'] < 0.05:
                                Q = 0.0
                                Sice = loaded_variables['Sice']
                                Sbr_ice = -Tice / mu
                                Vbr_ice = Sice / Sbr_ice
                                
                                Si = loaded_variables['Si']
                                Sbr_i = -(Ti - 273.15) / mu
                                Vbr_i = Si / Sbr_i
                                    
                                Sice_5 = loaded_variables['Sice_5']
                                Sbr_5 = Sice_5 / 0.05
                                Tice_5 = -mu * Sbr_5
                                    
                                hi_5 = (hi * loaded_variables['ki_5_bott'] * (Tice_5 - (Tfr - 273.15))) / ((loaded_variables['ki_ice_5'] * (Tice - Tice_5)) + (loaded_variables['ki_5_bott'] * (Tice_5 - (Tfr - 273.15))))
                                Tice_bio = (loaded_variables['ki_bio_bott'] * (Tfr - 273.15) + loaded_variables['ki_5_bio'] * Tice_5) / (loaded_variables['ki_bio_bott'] + loaded_variables['ki_5_bio'])
                                Sice_bio = loaded_variables['Sice_bio']
                                Sbr_bio = -Tice_bio / mu
                                Vbr_bio = Sice_bio / Sbr_bio
                                hi_bio = hi_5
                                    
                                Sice_bott = loaded_variables['Sice_bott']
                                Sbr_bott = -(Tfr - 273.1499) / mu
                                Vbr_bott = Sice_bott / Sbr_bott
                        
                        elif hi <= ZERO:
                            Q = ZERO
                            Sice = 0.0
                            Sice_bott = 0.0
                            Sbr_ice = 0.0
                            Vbr_ice = 0.0
                            
                            Si = 0.0
                            Sbr_i = 0.0
                            Vbr_i = 0.0
                            
                            Sbr_5 = 0.0
                            Sice_5 = 0.0
                            Tice_5 = Tfr - 273.15
                            hi_5 = ZERO
                            
                            Sbr_bio = 0.0
                            Sice_bio = 0.0
                            Tice_bio = Tfr - 273.15
                            Vbr_bio = 0.0
                            hi_bio = hi_5
                            
                            Sice_bott = 0.0
                            Sbr_bott = 0.0
                            Vbr_bott = 0.0
                            
                    elif hi <= ZERO:
                        strIII = '4'
                        ISI = Fs * ((infra_fr * np.exp(-k_ocean_red * h_mix)) + ((1 - infra_fr) * np.exp(-k_ocean_vis * h_mix)))
                        I0 = ISI
                       
                        Fsurf = 0.0
                        Fbott = 0.0
                            
                        IM = ISI
                        ISI_10 = ISI
                            
                        deltahs_melt_surf = 0.0
                        deltahmi_melt_surf = 0.0
                        deltahi_melt_surf = 0.0
                        deltahi_bott = 0.0
                        delta_bucket = 0.0
                        hs_prec_bucket = 0.0
                       
                        hs = ZERO
                        hi = ZERO
                        hmi = ZERO
                            
                        Ri = 0.0
                       
                        Q = ZERO
                        Sice = 0.0
                        Sice_bott = 0.0
                        Sbr_ice = 0.0
                        Vbr_ice = 0.0
                        Ssnowice = 0.0
                        
                        Si = 0.0
                        Sbr_i = 0.0
                        Vbr_i = 0.0
                       
                        Sbr_5 = 0.0
                        Sice_5 = 0.0
                        Tice_5 = Tfr - 273.15
                        hi_5 = ZERO
                       
                        Sbr_bio = 0
                        Sice_bio = 0.0
                        Tice_bio = Tfr - 273.15
                        Vbr_bio = 0.0
                        hi_bio = hi_5
                       
                        Sice_bott = 0.0
                        Sbr_bott = 0.0
                        Vbr_bott = 0.0
        else:
            print('ERROR MELTING SEASON')
            
        #% Properties depending on T and S
        Tice_av = Ti - 273.15
        Sice_av = (Sice_bott + Sice) / 2
            
        if Tice_av > -2.000:
            F1 = -4.1221e-2 - 1.8407 * 10 * Tice_av + 5.8402e-1 * Tice_av**2 + 2.1454e-1 * Tice_av**3
            F2 = 9.0312e-2 - 1.6111e-2 * Tice_av + 1.2291e-4 * Tice_av**2 + 1.3603e-4 * Tice_av**3
        else:
            F1 = -4.732 - 2.245 * 10 * Tice_av - 6.397e-1 * Tice_av**2 - 1.074e-2 * Tice_av**3
            F2 = 8.903e-2 - 1.763e-2 * Tice_av - 5.330e-4 * Tice_av**2 - 8.801e-6 * Tice_av**3
                
        roi_surf = 0.917 - 1.403e-4 * Tice  
        ro_sice_surf = 0.905  
        ro_sice_bott = 0.885  
        roi_av = 0.917 - 1.403e-4 * Tice_av
        ro_sice_bulk = (1 - Va) * ((roi_av * F1) / (F1 - (roi_av * Sice_av * F2)))

        ro_br = 1 + 8e-4 * Sbr_ice
        ro_br_i = 1 + 8e-4 * Sbr_i
        ro_br_5 = 1 + 8e-4 * Sbr_5
        ro_br_bio = 1 + 8e-4 * Sbr_bio
        ro_br_bott = 1 + 8e-4 * Sbr_bott

        k0i = 418.6 * (5.35e-3 - 2.568e-5 * Tice_av) 
        kb = 418.6 * (1.25e-3 + 3.0e-5 * Tice_av)

        k0i_5_bio = 418.6 * (5.35e-3 - 2.568e-5 * ((Tice_5 + Tice_bio) / 2))
        kb_5_bio = 418.6 * (1.25e-3 + 3.0e-5 * ((Tice_5 + Tice_bio) / 2))
        ki_5_bio = (1 - Va - 0.05) * k0i_5_bio + 0.005 * kb_5_bio

        k0i_ice_5 = 418.6 * (5.35e-3 - 2.568e-5 * ((Tice + Tice_5) / 2))
        kb_ice_5 = 418.6 * (1.25e-3 + 3.0e-5 * ((Tice + Tice_5) / 2))
        ki_ice_5 = (1 - Va - Vbr_ice) * k0i_ice_5 + Vbr_ice * kb_ice_5
            
        k0i_ice_bott = 418.6 * (5.35e-3 - 2.568e-5 * ((Tice + (Tfr - 273.15)) / 2))
        kb_ice_bott = 418.6 * (1.25e-3 + 3.0e-5 * ((Tice + (Tfr - 273.15)) / 2))
        ki_ice_bott = (1 - Va - Vbr_ice) * k0i_ice_bott + Vbr_ice * kb_ice_bott

        k0i_5_bott = 418.6 * (5.35e-3 - 2.568e-5 * ((Tice_5 + (Tfr - 273.15)) / 2))
        kb_5_bott = 418.6 * (1.25e-3 + 3.0e-5 * ((Tice_5 + (Tfr - 273.15)) / 2))
        ki_5_bott = (1 - Va - 0.05) * k0i_5_bott + 0.05 * kb_5_bott

        k0i_ice_bio = 418.6 * (5.35e-3 - 2.568e-5 * ((Tice + Tice_bio) / 2))
        kb_ice_bio = 418.6 * (1.25e-3 + 3.0e-5 * ((Tice + Tice_bio) / 2))
        ki_ice_bio = (1 - Va - Vbr_ice) * k0i_ice_bio + Vbr_ice * kb_ice_bio

        k0i_bio_bott = 418.6 * (5.35e-3 - 2.568e-5 * ((Tice_bio + (Tfr - 273.15)) / 2))
        kb_bio_bott = 418.6 * (1.25e-3 + 3.0e-5 * ((Tice_bio + (Tfr - 273.15)) / 2))
        ki_bio_bott = (1 - Va - Vbr_bio) * k0i_bio_bott + Vbr_bio * kb_bio_bott
            
        if hi > hi_min:
            if hi >= 0.1:
                ISI_layer = I0 * np.exp(-ksi_av * ((hi - 0.1) / 2))
                ISI_bio = I0 * np.exp(-ksi_av * (hi - 0.1 - hi_5 / 2))
            elif hi < 0.1:
                ISI_layer = IM * np.exp(-ksi_10_av * (hi / 2))
                ISI_bio = IM * np.exp(-ksi_10_av * (hi - hi_5 / 2))
        else:
            ISI_layer = 0.0
            ISI_bio = 0.0
                
        count_time_steps += 1

        index = loaded_variables['index'] + 1

        print(f'outcome of h1 in category {i}: ', "{:.6f}".format(hi))
        
        loaded_variables = {
            'index': index, 'Tmix': Tmix, 'T0': T0, 'Ts': Ts, 'Tsi': Tsi, 'Ti': Ti, 'hs': hs, 'hi': hi, 'hmi': hmi, 'hmi_new': hmi_new, 
            'hs_prec_bucket': hs_prec_bucket, 'hs_prec': hs_prec, 'hs_tot': hs_tot, 'delta_sublim': delta_sublim, 'delta_bucket': delta_bucket, 
            'snow_fr': snow_fr, 'sea_water_fr': sea_water_fr, 'ks': ks, 'ksi_10_av': ksi_10_av, 'ks_av': ks_av, 'ros': ros, 'ros_new': ros_new, 
            'Vbr_ice': Vbr_ice, 'Sbr_ice': Sbr_ice, 'Sbr_bott': Sbr_bott, 'Sice': Sice, 'Tice': Tice, 'Sice_bott': Sice_bott, 'Ssnowice': Ssnowice, 
            'Vbr_i': Vbr_i, 'Sbr_i': Sbr_i, 'Si': Si, 'Sice_5': Sice_5, 'Sbr_5': Sbr_5, 'Tice_5': Tice_5, 'hi_5': hi_5, 'Sice_bio': Sice_bio, 
            'Tice_bio': Tice_bio, 'Sbr_bio': Sbr_bio, 'hi_bio': hi_bio, 'Vbr_bio': Vbr_bio, 'Vbr_bott': Vbr_bott, 'ISI_bio': ISI_bio, 'R': R, 
            'ro_sice_bott': ro_sice_bott, 'ro_sice_surf': ro_sice_surf, 'ro_sice_bulk': ro_sice_bulk, 'ks_snow': ks_snow, 'alpha': alpha, 
            'ki_5_bott': ki_5_bott, 'ki_ice_5': ki_ice_5, 'ki_ice_bio': ki_ice_bio, 'ki_bio_bott': ki_bio_bott, 'ki_5_bio': ki_5_bio, 
            'ki_ice_bott': ki_ice_bott, 'ro_br': ro_br, 'ro_br_i': ro_br_i, 'ro_br_5': ro_br_5, 'ro_br_bio': ro_br_bio, 'ro_br_bott': ro_br_bott, 'F': F, 'ros_av': ros_av
        }

        globals().update(loaded_variables)

        # Save the variable to a new .pkl file with the same name
        output_filename = f'prev_index_values_{i}.pkl'

        # Print the name of the file that will be saved
        print(f"Saving file: {output_filename}")

        with open(output_filename, 'wb') as f:
            pickle.dump(loaded_variables, f)

        # Read the existing CSV file into a DataFrame
        df = pd.read_csv(csv_filepath)

        #for category in categories:
        column_name = f'h1_cat{i}'
        df.at[0, column_name] = hi

        #for category in categories:
        column_name = f'Fnet_cat{i}'
        df.at[0, column_name] = F
        
        # Save the DataFrame to a CSV file
        df.to_csv(csv_filepath, index=False)
