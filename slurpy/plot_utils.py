#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:22:51 2020

@author: jennywong
"""
import numpy as np
import matplotlib.pyplot as plt 

from slurpy.data_utils import readdata
from slurpy.getparameters import getcsbmassoxygen
from slurpy.coreproperties import icb_radius, aO
from slurpy.lookup import premdensity, liquidus

def plot_profile(inputDir):
    
# Read data
    try:
        print(inputDir)
        inputs,outputs,profiles = readdata(inputDir)
    except:
        print('No solution')
    radius=(profiles.r)*1e-3
    
    print('State = {}'.format(outputs.state.iloc[0]))
    
    # Plots
    w, h = plt.figaspect(0.75)*2
    # Temperature
    fig1=plt.figure
    fig1, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,figsize=(w,h))
    ax1.plot(radius,profiles.temp,color='orange')
    ax1.set(ylabel="Temperature (K)")
    
    # Oxygen
    (mass_conc_O,acore) =getcsbmassoxygen(inputs.oxygen_bulk)
    acore=float(acore)
    ax2.plot(radius,profiles.oxygen*acore/aO*100,color='orange')
    ax2.set(ylabel="Oxygen (mol.%)")
    
    # Solid flux
    ax3.plot(radius,profiles.solidflux,color='orange')
    ax3.set(xlabel="Radius (km)",ylabel="Solid flux ($\mathrm{kg m^{-2} s^{-1}}$)")
    
    # Density
    ax4.plot(radius,profiles.density,color='orange')
    # PREM   
    density_prem=premdensity(radius*1e3)
    ax4.plot(radius,density_prem,'k--')
    ax4.legend(["Spherical", "PREM"])
    
    ax4.legend(["Spherical", "PREM"])
    ax4.set(xlabel="Radius (km)",ylabel="Density ($\mathrm{kg m^{-3}}$)")
    
    # Liquidus        
    radius_liquidus=np.linspace(icb_radius,(icb_radius+400e3))
    temp_liquidus=liquidus(radius_liquidus)
    ax1.plot(radius_liquidus*1e-3,temp_liquidus,'k--')
    
    # Constant CSB oxygen
    ax2.hlines(y=8, xmin=radius[0], xmax=icb_radius*1e-3+400,colors='k',linestyles='dashed')
           
    # Zero solid flux
    ax3.hlines(y=0, xmin=radius[0], xmax=icb_radius*1e-3+400,colors='k',linestyles='dashed')
     
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
