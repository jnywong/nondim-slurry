#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 14:52:41 2020

@author: jennywong

###############################################################################
# script_sensitivity.py                                                       #
###############################################################################

# Perform sensitivity study outlined in Wong et al. (in prep).

Parameters
----------
plotOn : int
    Toggle 0='off' and 1='on' to display the solution output. Recommended off for large parameter searches but useful for probing a small number of solutions
layer_thicknesses : numpy array
    Specify slurry layer thickness (m).
thermal_conductivities : numpy array
    Specify thermal conductivity (W/m/K).
icb_heatfluxes : numpy array
    Specify heat flux through the ICB (TW).
csb_heatfluxes : numpy array
    Specify heat flux through the CSB (TW).
h : float
    Specify step size through the heat fluxes for parameter searches.
sensitivityOn : int
    Toggle manual control of the CSB temperature and/or CSB oxygen concentration for sensitivity studies.
mol_conc_oxygen_bulk : float
    If sensitivityOn=1, then manually specify the CSB oxygen concentration. Default value is 8 mol.%.
csb_temp : float
    If sensitivityOn=1, then manually specify the CSB temperture (K). Default value depends on the layer thickness, given by function slurpy.lookup.liquidus.
sedimentation_constant : float
    If sensitivityOn=1, then manually specify the sedimentation constant used to derive the solid fraction from the solid flux. Default value is 1e-2 kg s/m^3.

Returns
-------
outputDir : str
    Folder name solution outputs are saved in
radius : numpy array
    Solution mesh (m)
temp : numpy array
    Solution temperature (K)
xi : numpy array
    Solution oxygen concentration (mass and not molar concentration)
solidFlux : numpy array
    Solution solid flux (kg/m^2/s)
density : numpy array
    Solution density (kg/m^3)

"""

# %% IMPORT STATEMENTS

import numpy as np
import pickle
import os

from slurpy.slurry import solveslurry
from slurpy.lookup import liquidus
from slurpy.getparameters import getcsbradius
from slurpy.plot_utils import plot_sensitivity, plot_sedimentation

# %% MODEL INPUTS
# Save plot?
saveOn=0

# Input parameters
layer_thickness=150e3 # (m)
thermal_conductivity=100. # (W m^-1 K^-1)
icb_heatflux=2.5 # (TW)
csb_heatflux=5.0 # (TW)
h=0.05 # stepsize of heat flux through parameter space

# Sensitivity study
csb_temp = np.arange(4500.,6100.,100) # (K)
csb_oxy = np.arange(2,12.5,0.5) # (mol.%)
sed_con= np.array([1e-5,1e-4,1e-3,1e-2,1e-1]) # (kg s/m^3) pre-factor in sedimentation coefficient, b(phi)
#------------------------------------------------------------------------------
# %% RUN THE CODE
# Append default values
csb_radius=getcsbradius(layer_thickness)
csb_temp0=liquidus(csb_radius)
csb_temp=np.append(csb_temp,csb_temp0)
csb_temp=np.sort(csb_temp)
csb_oxy0=8.
if csb_oxy0 not in csb_oxy:
    csb_oxy=np.append(csb_oxy,csb_oxy0)
    csb_oxy=np.sort(csb_oxy)

n_temp = csb_temp.size
n_oxy = csb_oxy.size
n_sed = sed_con.size

k=1 # counter to track progress
# Make sensitivity directory
if not os.path.exists('sensitivity'):
    os.makedirs('sensitivity')

# Loop through csb temperature
for x in [x for x in csb_temp]:
    filename = 'sensitivity/temp_{:.0f}'.format(x).replace('.','_')
    if os.path.exists(filename):
        continue # skip if file exists
    (outputDir,radius,temp,xi,solidFlux,density)= \
        solveslurry(layer_thickness, icb_heatflux, csb_heatflux, \
                    thermal_conductivity, x, h)
    print('CSB temp: Run {}/{}'.format(k,n_temp))

    with open(filename,'wb') as f:
        pickle.dump([radius,temp,xi,solidFlux,density], f)
    print('{} saved'.format(filename))

    k+=1

k=1
for y in [y for y in csb_oxy]:
    filename = 'sensitivity/oxy_{:.1f}'.format(y).replace('.','_')
    if os.path.exists(filename):
        continue # skip if file exists
    (outputDir,radius,temp,xi,solidFlux,density)= \
        solveslurry(layer_thickness, icb_heatflux, csb_heatflux, \
                    thermal_conductivity, csb_temp0, h, mol_conc_oxygen_bulk=y)
    print('CSB oxy: Run {}/{}'.format(k,n_oxy))

    with open(filename,'wb') as f:
        pickle.dump([radius,temp,xi,solidFlux,density], f)
    print('{} saved'.format(filename))

    k+=1

plot_sensitivity(csb_temp,csb_oxy,csb_temp0,csb_oxy0,saveOn)

k=1
for z in [z for z in sed_con]:
    filename = 'sensitivity/sed_{:.0f}'.format(np.log10(z)).replace('.','_')
    if os.path.exists(filename):
        continue # skip if file exists
    (outputDir,radius,temp,xi,solidFlux,density)= \
        solveslurry(layer_thickness, icb_heatflux, csb_heatflux, \
                    thermal_conductivity, csb_temp0, h, sedimentation_constant=z)
    print('Sedimentation constant: Run {}/{}'.format(k,n_sed))

    with open(filename,'wb') as f:
        pickle.dump([radius,temp,xi,solidFlux,density], f)
    print('{} saved'.format(filename))

    k+=1

plot_sedimentation(sed_con,saveOn,mol_conc_oxygen_bulk=8,figAspect=0.75)
