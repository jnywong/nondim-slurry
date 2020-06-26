#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:20:26 2019

@author: Jenny Wong

###############################################################################
# script_parametersearch.py                                                   #
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

from slurpy.slurry import solveslurry
from slurpy.getparameters import getcsbradius
from slurpy.lookup import liquidus
from slurpy.plot_utils import plot_profile

# %% MODEL INPUTS
# Show plots?
plotOn=1 # show temp, xi, solid flux and density profiles

overwriteOn=1

# Input parameters
layer_thicknesses=np.array([250e3]) # (m)
# layer_thicknesses=np.array([150e3,200e3,250e3,300e3,350e3,400e3]) #(m)
thermal_conductivities=np.array([100.]) # (W m^-1 K^-1)
icb_heatfluxes=np.array([2.05]) # (TW)
csb_heatfluxes=np.array([6.55]) # (TW)

h=0.05 # stepsize of heat flux through parameter space
# csb_heatfluxes=np.arange(5,6.05,h) # (TW)
# icb_heatfluxes=np.arange(2,4.05,h) # (TW)
# csb_heatfluxes=np.arange(0.05,15.05,h) # (TW)
# icb_heatfluxes=np.arange(0.05,5.05,h) # (TW)

#------------------------------------------------------------------------------
# %% RUN THE CODE
n_thick = layer_thicknesses.size
n_thermal = thermal_conductivities.size
n_icb = icb_heatfluxes.size
n_csb = csb_heatfluxes.size
n_tot = n_thick*n_thermal*n_icb*n_csb
k=1 # counter to track parameter search progress

for w,x,y,z in [(w,x,y,z) for w in layer_thicknesses for x in icb_heatfluxes for y in csb_heatfluxes for z in thermal_conductivities]:
    print('Run {}/{}: d = {:.0f}, Qi = {}, Qsl = {}, k = {:.0f}'.format(k,n_tot,w*1e-3,x,y,z))
    csb_radius = getcsbradius(w)
    csb_temp = liquidus(csb_radius)
    (outputDir,radius,temp,xi,solidFlux,density)=solveslurry(w,x,y,z,csb_temp,h,overwrite=overwriteOn)
    k+=1

    # %%PLOT
    if plotOn==1:
        plot_profile(outputDir)
