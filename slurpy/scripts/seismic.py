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
perseuscsb_heatflux : numpy float
    Specify heat flux through the CSB (TW).

Returns
-------
foldername : str
    Folder name of Lewis number solution is saved in
filename : str
    Filename solution is saved in

"""

# %% IMPORT STATEMENTS
import numpy as np
from slurpy.plot_utils import plot_seismic
from slurpy.coreproperties import density_solidFe

# %% MODEL INPUTS
# Save plot?
saveOn=1

# Input parameters
layer_thickness=np.array([150e3]) # (m)
thermal_conductivity=np.array([100.]) # (W m^-1 K^-1)
icb_heatflux=np.array([2.5]) # (TW)
csb_heatflux=np.array([5])
seis='prem'

#------------------------------------------------------------------------------
# %% RUN THE CODE
# Load solution
# foldername, filename = get_outputDir(layer_thickness,icb_heatflux,csb_heatflux, \
#                                      thermal_conductivity, model=seis)
# Plot and save
radius,slurry_vp,slurry_density,ohtaki_vp = plot_seismic(layer_thickness, thermal_conductivity,
                                          icb_heatflux, csb_heatflux, saveOn)
