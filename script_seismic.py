#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 1 14:52:41 2020

@author: jennywong

###############################################################################
# script_seismic.py                                                           #
###############################################################################

# Compare slurry P wave speed with PREM, ak135 and Ohtaki et al. (2015) as 
outlined in Wong et al. (in prep).

Parameters
----------
saveOn : int
    Toggle 0='off' and 1='on' to save figure output.
layer_thickness : numpy float
    Specify slurry layer thickness (m).
thermal_conductivity : numpy float
    Specify thermal conductivity (W/m/K).
icb_heatflux : numpy float
    Specify heat flux through the ICB (TW).
csb_heatflux : numpy float
    Specify heat flux through the CSB (TW).

Returns
-------
foldername : str
    Folder name of Lewis number solution is saved in
filename : str
    Filename solution is saved in

"""

# %% IMPORT STATEMENTS
from slurpy.data_utils import get_outputDir
from slurpy.plot_utils import plot_seismic

# %% MODEL INPUTS
# Save plot?
saveOn=1

# Input parameters
layer_thickness=250e3 # (m)
thermal_conductivity=100. # (W m^-1 K^-1)
icb_heatflux=3. # (TW)
csb_heatflux=6. # (TW)
#------------------------------------------------------------------------------
# %% RUN THE CODE
# Load solution
foldername, filename = get_outputDir(layer_thickness, icb_heatflux, \
                                     csb_heatflux, thermal_conductivity)
# Plot and save
plot_seismic(foldername, filename, saveOn)
