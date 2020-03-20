#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:20:26 2019

@author: Jenny Wong

"""

###############################################################################
# NONDIMENSIONAL SPHERICAL STEADY STATE SLURRY                                #
###############################################################################
#
# This script allows the user to solve the nondimensional spherical steady
# state slurry equations.
#

# %% IMPORT STATEMENTS

import numpy as np
import pickle
import os

# from mainroutines import solveslurry
from slurpy import slurry
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

# Default parameters
sensitivityOn=0 # switch for sensitivity study
mol_conc_oxygen_bulk = 8
csb_temp = 5000
sedimentation_constant=1e-2 # pre-factor in sedimentation coefficient, b(phi)
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
    (outputDir,radius,temp,xi,solidFlux,density)= slurry.solveslurry(w,x,y,z, \
                       mol_conc_oxygen_bulk,csb_temp,sedimentation_constant)
    
    print('Run {}/{}'.format(k,n_tot))
    k=k+1

    # SENSITIVITY
    if sensitivityOn==1:
        filename = 'sensitivity/temp_{:.0f}'.format(csb_temp).replace('.','_')
        filename = 'sensitivity/xi_{:.1f}'.format(mol_conc_oxygen_bulk).replace('.','_')   
        if not os.path.exists(filename):
            os.makedirs(filename)
        with open(filename,'wb') as f:
            pickle.dump([radius,temp,xi,solidFlux,density], f)
    
    # PLOT
    if plotOn==1:
        plot_profile(outputDir)