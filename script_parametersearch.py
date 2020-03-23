#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:20:26 2019

@author: Jenny Wong

###############################################################################
# NONDIMENSIONAL SPHERICAL STEADY STATE SLURRY                                #
###############################################################################

# Solve the 1D, steady, spherical slurry system outlined in Wong et al. (in prep)
# (see also Wong et al. 2018).

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

# from mainroutines import solveslurry
from slurpy.slurry import solveslurry
from slurpy.lookup import liquidus
from slurpy.getparameters import getcsbradius
from slurpy.plot_utils import plot_profile

# %% MODEL INPUTS
# Show plots?
plotOn=1 # show temp, xi, solid flux and density profiles

# Input parameters
layer_thicknesses=np.array([150e3]) # (m)
#layer_thicknesses=np.array([150e3,200e3,250e3,300e3,350e3,400e3]) #(m)
thermal_conductivities=np.array([100.0]) # (W m^-1 K^-1)
icb_heatfluxes=np.array([2.5]) # (TW)
csb_heatfluxes=np.array([5.5]) # (TW)

h=0.05 # stepsize of heat flux through parameter space
# csb_heatfluxes=np.arange(0.05,15.05,h) # (TW)
# icb_heatfluxes=np.arange(0.05,5.05,h) # (TW)

# Sensitivity study
sensitivityOn=0 # switch for sensitivity study
mol_conc_oxygen_bulk = 8 #(mol.%)
csb_temp = 5000 #(K)
sedimentation_constant=1e-2 # (kg s/m^3) pre-factor in sedimentation coefficient, b(phi)
#------------------------------------------------------------------------------
# %% RUN THE CODE
n_icb = icb_heatfluxes.size
n_csb = csb_heatfluxes.size
n_tot = n_icb*n_csb
k=1 # counter to track parameter search progress

for w,x,y,z in [(w,x,y,z) for w in layer_thicknesses for x in icb_heatfluxes for y in csb_heatfluxes for z in thermal_conductivities]:
    if sensitivityOn ==0:
        csb_radius = getcsbradius(w)
        csb_temp = liquidus(csb_radius)
    elif sensitivityOn ==1:
        filename_temp = 'sensitivity/temp_{:.0f}'.format(csb_temp).replace('.','_')
        filename_oxy = 'sensitivity/xi_{:.1f}'.format(mol_conc_oxygen_bulk).replace('.','_')
        if not os.path.exists(filename_temp):
            os.makedirs(filename_temp)
        if not os.path.exists(filename_oxy):
            os.makedirs(filename_oxy)            
    (outputDir,radius,temp,xi,solidFlux,density,csb_temp)= solveslurry(w,x,y,z, \
                       mol_conc_oxygen_bulk,csb_temp,sedimentation_constant)
    if sensitivityOn ==1:
        with open(filename_temp,'wb') as f:
            pickle.dump([radius,temp,xi,solidFlux,density], f)
        with open(filename_oxy,'wb') as f:
            pickle.dump([radius,temp,xi,solidFlux,density], f)

    print('Run {}/{}'.format(k,n_tot))
    k=k+1

    # PLOT
    if plotOn==1:
        plot_profile(outputDir)
