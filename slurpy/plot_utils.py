#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:22:51 2020

@author: jennywong
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import pickle
import os

from slurpy.data_utils import readdata
from slurpy.getparameters import getcsbmassoxygen
from slurpy.coreproperties import icb_radius, aO
from slurpy.lookup import premdensity, liquidus

# %% 
def plot_profile(inputDir):
    
# Read data
    try:
        # print(inputDir)
        inputs,outputs,profiles = readdata(inputDir)
    except:
        print('No solution')
    
    # print('State = {}'.format(outputs.state.iloc[0]))
    
    # Plots
    csb_radius = pd.to_numeric(profiles.r.iloc[-1])
    radius=(profiles.r)*1e-3
    w, h = plt.figaspect(0.75)*2
    # Temperature
    fig1=plt.figure
    fig1, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,figsize=(w,h))
    ax1.plot(radius,profiles.temp)
    ax1.set(ylabel="Temperature (K)")
    
    # Oxygen
    (mass_conc_O,acore) =getcsbmassoxygen(inputs.oxygen_bulk)
    acore=float(acore)
    ax2.plot(radius,profiles.oxygen*acore/aO*100)
    ax2.set(ylabel="Oxygen (mol.%)")
    
    # Solid flux
    ax3.plot(radius,profiles.solidflux)
    ax3.set(xlabel="Radius (km)",ylabel="Solid flux ($\mathrm{kg m^{-2} s^{-1}}$)")
    
    # Density
    ax4.plot(radius,profiles.density)
    # PREM   
    density_prem=premdensity(radius*1e3)
    ax4.plot(radius,density_prem,'k--',label='PREM')
    ax4.legend()

    ax4.set(xlabel="Radius (km)",ylabel="Density ($\mathrm{kg m^{-3}}$)")
    
    # Liquidus        
    radius_liquidus=np.linspace(icb_radius,csb_radius)
    temp_liquidus=liquidus(radius_liquidus)
    ax1.plot(radius_liquidus*1e-3,temp_liquidus,'k--', label = 'Liquidus (Davies et al. 2015)')
    ax1.legend()
    
    # Constant CSB oxygen
    # ax2.hlines(y=8, xmin=radius[0], xmax=radius.iloc[-1],colors='k',linestyles='dashed')
           
    # Zero solid flux
    # ax3.hlines(y=0, xmin=radius[0], xmax=radius.iloc[-1],colors='k',linestyles='dashed')
     
    # PREM
    radius_prem=np.linspace(icb_radius,csb_radius)
    density_prem=premdensity(radius_prem)
    ax4.plot(radius_prem*1e-3,density_prem,'k--')
    ax4.set(xlabel="Radius (km)",ylabel="Density ($\mathrm{kg m^{-3}}$)")
    
    ax1.set(ylabel="Temperature (K)")
    ax2.set(ylabel="Oxygen (mol.%)")
    ax3.set(xlabel="Radius (km)",ylabel="Solid flux ($\mathrm{kg m^{-2} s^{-1}}$)")
    
    ax1.set_xlim([radius.iloc[0],radius.iloc[-1]])
    
    plt.tight_layout()
    
    plt.show()

# %%
def plot_sensitivity(csb_temp,csb_oxygen,csb_temp0,csb_oxy0,saveOn,aspectRatio=0.75):
    w, h = plt.figaspect(aspectRatio)
    fig1, (ax1,ax2) = plt.subplots(1,2,figsize=(2*w,h),sharey=True)
    
    # Temperature
    nTemp = csb_temp.size
    colors=plt.cm.OrRd(np.linspace(0.4,1,nTemp))
    den_jump=[]
    
    for i in range(nTemp):
        filename = 'sensitivity/temp_{:.0f}'.format(csb_temp[i]).replace('.','_')
        with open(filename, 'rb') as f:
            (radius,temp,xi,solidFlux,density)=pickle.load(f)
        if i ==0 or i == nTemp-1:
            ax1.plot(radius*1e-3,density,color=colors[i], linewidth = 2,
                     label =r'$T_l=${:.0f} K'.format(csb_temp[i]))
        # Reference case
        elif csb_temp[i]==csb_temp0:
            den_jump0 = density[0]-density[-1]
            ax1.plot(radius*1e-3,density,color='silver', linewidth = 2,
                     label=r'$T_l=$5457 K')
        else:
            ax1.plot(radius*1e-3,density,color=colors[i])
        den_jump.append(density[0]-density[-1])
    
    # PREM
    density_prem=premdensity(radius)
    ax1.plot(radius*1e-3,density_prem, 'k', linestyle = '--', label = r'PREM')
    
    ax1.legend(fontsize=11.5)
    ax1.set_xlim([radius[0]*1e-3,radius[-1]*1e-3])
    ax1.set(xlabel="Radius (km)",ylabel="Density ($\mathrm{kg m^{-3}}$)")
    #ax1.set_xticks(np.arange(1225,1380,25))
    
    den_low = (den_jump[0]-den_jump0)/den_jump0*100
    den_high = (den_jump[-1]-den_jump0)/den_jump0*100
    print('Temperature: Density jump ranges from {:.2f}% to {:.2f}% of reference'.format(den_low, den_high))
    
    # Oxygen
    nOxy = csb_oxygen.size
    colors=plt.cm.GnBu(np.linspace(0.4,1,nOxy))
    den_jump=[]
    
    for i in range(nOxy):
        filename = 'sensitivity/oxy_{:.1f}'.format(csb_oxygen[i]).replace('.','_')
        with open(filename, 'rb') as f:
            (radius,temp,xi,solidFlux,density)=pickle.load(f)
        if i ==0 or i == nOxy-1:
            ax2.plot(radius*1e-3,density,color=colors[i], linewidth = 2,
                     label =r'$\xi_O=${:.1f} mol.%'.format(csb_oxygen[i]))
        # Reference case
        elif csb_oxygen[i]==csb_oxy0:
            den_jump0 = density[0]-density[-1]
            ax2.plot(radius*1e-3,density,color='silver', linewidth = 2,
                     label=r'$\xi_O=$8.0 mol.%')
        else:
            ax2.plot(radius*1e-3,density,color=colors[i])
        den_jump.append(density[0]-density[-1])
    
    # PREM
    density_prem=premdensity(radius)
    ax2.plot(radius*1e-3,density_prem, 'k', linestyle = '--')
    
    ax2.legend(fontsize=11.5)
    ax2.set_xlim([radius[0]*1e-3,radius[-1]*1e-3])
    # ax2.set_ylim([12080,max_den])
    ax2.set(xlabel="Radius (km)")
    
    den_low = (den_jump[0]-den_jump0)/den_jump0*100
    den_high = (den_jump[-1]-den_jump0)/den_jump0*100
    print('Oxygen: Density jump ranges from {:.2f}% to {:.2f}% of reference'.format(den_low, den_high))
    
    if saveOn==1:
        saveDir='figures/sensitivity/'
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        fig1.savefig(saveDir+"temp_oxy.pdf",format='pdf', dpi=200, bbox_inches='tight')
        fig1.savefig(saveDir+"temp_oxy.png",format='png', dpi=200, bbox_inches='tight')
        print('Figure saved as {}'.format(saveDir+"temp_oxy.pdf"))
