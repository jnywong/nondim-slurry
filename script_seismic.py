#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 1 14:52:41 2020

@author: jennywong

###############################################################################
# script_vp.py                                                       #
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

from slurpy.data_utils import get_outputDir, readdata

from slurpy.slurry import solveslurry
from slurpy.lookup import vpspeed,premdensity
from slurpy.getparameters import getcsbradius
from slurpy.plot_utils import plot_seismic

# %% MODEL INPUTS
# Save plot?
saveOn=1 

# Input parameters
layer_thickness=250e3 # (m)
thermal_conductivity=100. # (W m^-1 K^-1)
icb_heatflux=3.5 # (TW)
csb_heatflux=6. # (TW)
#------------------------------------------------------------------------------
# %% RUN THE CODE

# Make vp directory
if not os.path.exists('seismic'):
    os.makedirs('seismic')

# Load solution
inputDir = get_outputDir(layer_thickness, icb_heatflux, csb_heatflux, thermal_conductivity)

# Plot
plot_seismic(inputDir,saveOn)
