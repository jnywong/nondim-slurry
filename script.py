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
# _initialise: use previous solution to initialise new solution
#
# This script allows the user to solve the nondimensional spherical steady
# state slurry equations outlined in the document "A nondimensional steady
# slurry model".
#
# Parameters
# ----------
# 
#
# Returns
# -------
#

# %% IMPORT STATEMENTS

import numpy as np
#matplotlib.use('PS')
import matplotlib.pyplot as plt
import pickle

from mainroutines import solveslurry
from readscripts import read_sphdata
from lookup import premdensity, liquidus
from coreproperties import aO, icb_radius
from refparameters import getcsbmassoxygen

# %% MODEL INPUTS

# Show plots?
plotOn=0 # plot figures

# Input parameters
layer_thicknesses=np.array([150e3]) #(km)
thermal_conductivities=np.array([100.0]) # (W m^-1 K^-1)
# icb_heatfluxes=np.array([2.5]) # (TW)
#csb_radius = getcsbradius(layer_thicknesses[0])
#csb_heatflux = icb_heatfluxes[0]*csb_radius**2/icb_radius**2
#csb_heatfluxes=np.array([csb_heatflux]) # (TW)
# csb_heatfluxes=np.array([5]) # (TW)

#layer_thicknesses=np.array([150e3,200e3,250e3,300e3,350e3,400e3]) #(m)

h=0.05 # stepsize of heat flux through parameter space
csb_heatfluxes=np.arange(0.05,15.05,h) # (TW) # low k
icb_heatfluxes=np.arange(0.05,5.05,h) # (TW)

# Default parameters
mol_conc_oxygen_bulk = 8
mol_conc_SSi = 8
self_diffusion= 0.98e-8
#core_cooling_rate= -100 # (K/Ga)
sedimentation_constant=1e-2 # pre-factor in sedimentation coefficient, b(phi)
initial_F=5. # initial guess for eigenvalue in boundary value problem
ic_age= 0.5 # inner core age (Ga), initial guess for snow speed
n=100 # mesh resolution
#------------------------------------------------------------------------------

# %% RUN THE CODE
n_icb = icb_heatfluxes.size
n_csb = csb_heatfluxes.size
n_tot = n_icb*n_csb
k=0
for w,x,y,z in [(w,x,y,z) for w in layer_thicknesses for x in icb_heatfluxes for y in csb_heatfluxes for z in thermal_conductivities]:    
    (outputDir,radius,temp,xi,solidFlux,density,csb_temp)=solveslurry(w,x,y,z,ic_age,  \
         self_diffusion,sedimentation_constant, \
         mol_conc_oxygen_bulk,mol_conc_SSi, \
         initial_F,n,h)
    
    print('Run {}/{}'.format(k,n_tot))
    k=k+1
    
#     ,radius,temp,xi,solidFlux,density,csb_temp)

#    # SENSITIVITY
#    filename = 'sensitivity/temp_{:.0f}'.format(csb_temp).replace('.','_')
#    filename = 'sensitivity/xi_{:.1f}'.format(mol_conc_oxygen_bulk).replace('.','_')    
#    with open(filename,'wb') as f:
#        pickle.dump([radius,temp,xi,solidFlux,density], f)
    
    if plotOn==1:
        # Read outputs
        try:
            print(outputDir)
            sph_inputs,sph_outputs,sph_profiles = read_sphdata(outputDir)
        except:
            print('No solution')
            continue
        sph_radius=(sph_profiles.z)*1e-3
        
        print('State = {}'.format(sph_outputs.state.iloc[0]))
        
        # Plots
        w, h = plt.figaspect(0.75)*2
        # Temperature
        fig1=plt.figure
        fig1, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,figsize=(w,h))
        ax1.plot(sph_radius,sph_profiles.temp,color='orange')
#        ax1.plot(radius_analytic*1e-3,temp_analytic,color='blue')
        ax1.set(ylabel="Temperature (K)")
        
        # Oxygen
        (mass_conc_O,acore) =getcsbmassoxygen(mol_conc_oxygen_bulk, mol_conc_SSi)
        acore=float(acore)
        aO=float(aO)
        ax2.plot(sph_radius,sph_profiles.oxygen*acore/aO*100,color='orange')
        ax2.set(ylabel="Oxygen (mol.%)")
        
        # Solid flux
        ax3.plot(sph_radius,sph_profiles.solidflux,color='orange')
        ax3.set(xlabel="Radius (km)",ylabel="Solid flux ($\mathrm{kg m^{-2} s^{-1}}$)")
        
        # Density
        ax4.plot(sph_radius,sph_profiles.density,color='orange')
        # PREM   
        density_prem=premdensity(sph_radius*1e3)
        ax4.plot(sph_radius,density_prem,'k--')
        ax4.legend(["Spherical", "PREM"])
        
        ax4.legend(["Spherical", "PREM"])
        ax4.set(xlabel="Radius (km)",ylabel="Density ($\mathrm{kg m^{-3}}$)")
        
        # Liquidus        
        radius_liquidus=np.linspace(icb_radius,(icb_radius+400e3))
        temp_liquidus=liquidus(radius_liquidus)
        ax1.plot(radius_liquidus*1e-3,temp_liquidus,'k--')
        
        # Constant CSB oxygen
        ax2.hlines(y=8, xmin=sph_radius[0], xmax=icb_radius*1e-3+400,colors='k',linestyles='dashed')
               
        # Zero solid flux
        ax3.hlines(y=0, xmin=sph_radius[0], xmax=icb_radius*1e-3+400,colors='k',linestyles='dashed')
         
        # PREM
        radius_prem=np.linspace(icb_radius,icb_radius+400e3)
        density_prem=premdensity(radius_prem)
        ax4.plot(radius_prem*1e-3,density_prem,'k--')
        ax4.set(xlabel="Radius (km)",ylabel="Density ($\mathrm{kg m^{-3}}$)")
        
        ax1.set(ylabel="Temperature (K)")
        ax2.set(ylabel="Oxygen (mol.%)")
        ax3.set(xlabel="Radius (km)",ylabel="Solid flux ($\mathrm{kg m^{-2} s^{-1}}$)")
        
        ax1.set_xlim([1220,icb_radius*1e-3+400])
        
        plt.tight_layout()
        
        plt.show()
#