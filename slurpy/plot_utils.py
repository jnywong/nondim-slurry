#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:22:51 2020

@author: jennywong
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pickle
import os
import matplotlib as mpl
mpl.rcParams['text.usetex'] = False

from slurpy.data_utils import readdata
from slurpy.getparameters import getcsbmassoxygen, getKphi, getphi
from slurpy.coreproperties import icb_radius, earth_radius, aO
from slurpy.lookup import premdensity, liquidus, premvp, ak135radius, ak135vp

# %% 
def plot_profile(inputDir):
    
# Read data
    try:
        # print(inputDir)
        inputs,outputs,profiles = readdata(inputDir)
    except:
        print('Folder does not exist')
    
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
                     label =r'$T_{{sl}}=${:.0f} K'.format(csb_temp[i]))
        # Reference case
        elif csb_temp[i]==csb_temp0:
            den_jump0 = density[0]-density[-1]
            ax1.plot(radius*1e-3,density,color='silver', linewidth = 2,
                     label=r'$T_\mathregular{{sl}}=$5457 K')
        else:
            ax1.plot(radius*1e-3,density,color=colors[i])
        den_jump.append(density[0]-density[-1])
    
    # PREM
    density_prem=premdensity(radius)
    ax1.plot(radius*1e-3,density_prem, 'k', linestyle = '--')
    
    ax1.legend(fontsize=11.5)
    ax1.set_xlim([radius[0]*1e-3,radius[-1]*1e-3])
    ax1.set(xlabel="Radius (km)",ylabel="Density ($\mathrm{kg m^{-3}}$)")
    
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
                     label =r'$\xi_{{sl}}=${:.1f} mol.%'.format(csb_oxygen[i]))
        # Reference case
        elif csb_oxygen[i]==csb_oxy0:
            den_jump0 = density[0]-density[-1]
            ax2.plot(radius*1e-3,density,color='silver', linewidth = 2,
                     label=r'$\xi_{{sl}}=$8.0 mol.%')
        else:
            ax2.plot(radius*1e-3,density,color=colors[i])
        den_jump.append(density[0]-density[-1])
    
    # PREM
    density_prem=premdensity(radius)
    ax2.plot(radius*1e-3,density_prem, 'k', linestyle = '--', label = r'PREM')
    
    ax2.legend(fontsize=11.5)
    ax2.set_xlim([radius[0]*1e-3,radius[-1]*1e-3])
    ax2.set(xlabel="Radius (km)")
    
    den_low = (den_jump[0]-den_jump0)/den_jump0*100
    den_high = (den_jump[-1]-den_jump0)/den_jump0*100
    print('Oxygen: Density jump ranges from {:.2f}% to {:.2f}% of reference'.format(den_low, den_high))
    
    ax1.set_title('(a)',x=0.95,y=1,fontsize=14)
    ax2.set_title('(b)',x=0.95,y=1,fontsize=14)    
    
    if saveOn==1:
        saveDir='figures/sensitivity/'
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        fig1.savefig(saveDir+"temp_oxy.pdf",format='pdf', dpi=200, bbox_inches='tight')
        fig1.savefig(saveDir+"temp_oxy.png",format='png', dpi=200, bbox_inches='tight')
        print('Figure saved as {}'.format(saveDir+"temp_oxy.pdf"))
        
# %%
def plot_sedimentation(sed_con,saveOn,mol_conc_oxygen_bulk=8,figAspect=0.75):
    nSed = sed_con.size
    colors=plt.cm.copper_r(np.linspace(0.4,1,nSed))
    
    w, h = plt.figaspect(figAspect)
    fig, (ax1,ax2) = plt.subplots(2,1,sharex=True,figsize=(w,2*h))
    
    # Format function for scientific notation in legend
    my_fun = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    
    for i in range(nSed):
        filename = 'sensitivity/sed_{:.0f}'.format(np.log10(sed_con[i])).replace('.','_')
        with open(filename, 'rb') as f:
            (radius,temp,xi,solidFlux,density)=pickle.load(f)
        # FIX:
        # ax1.plot(radius*1e-3,density,label=r'$k_\phi =$ {} $\mathrm{{kgm^{{-3}}s}}$'.format(my_fun.format_data(sed_con[i])),
        #        color=colors[i]) 
        ax1.plot(radius*1e-3,density,label=r'$k_\phi = {} \mathrm{{kgm^{{-3}}s}}$'.format(my_fun.format_data(sed_con[i])).replace('{','{{').replace('}','}}'),
               color=colors[i]) 
        ax1.plot(radius*1e-3,density,color=colors[i]) 
        Kphi = getKphi(sed_con[i],radius,mol_conc_oxygen_bulk)
        phi = getphi(Kphi,solidFlux)
        ax2.plot(radius*1e-3,phi,color=colors[i])
                 
    # PREM       
    density_prem=premdensity(radius)
    ax1.plot(radius*1e-3,density_prem,'k--', label='PREM')
    ax1.set(ylabel="Density ($\mathrm{kg m^{-3}}$)") #,yscale="log")
    ax2.set(xlabel="Radius (km)",ylabel="Solid fraction",yscale='log')
    ax2.axhline(0.6,color='k',linestyle='--') # rheological transition
    # labels = ['$k_\phi =${} $\mathrm{{kgm^{{-3}}s}}$'.format(my_fun.format_data(sed_con[i])),
              # '1','2','3','4']
    fig.legend(loc='center right', bbox_to_anchor=(1.4, 0.5),fontsize = 11.5)
    ax1.set_xlim([radius[0]*1e-3,radius[-1]*1e-3])
    ax2.set_ylim([1e-4,1])
    
    ax1.set_title('(a)',x=0.95,y=1,fontsize=14)
    ax2.set_title('(b)',x=0.95,y=1,fontsize=14)
    
    if saveOn==1:
        saveDir='figures/sensitivity/'
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        fig.savefig(saveDir+"sedimentation.pdf",format='pdf', dpi=200, bbox_inches='tight')
        fig.savefig(saveDir+"sedimentation.png",format='png', dpi=200, bbox_inches='tight')
        print('Figure saved as {}'.format(saveDir+"sedimentation.pdf"))
    plt.show()

# %%
def plot_seismic(foldername,filename,saveOn,figAspect=0.75):
    w, h = plt.figaspect(figAspect)
    fig, ax = plt.subplots(1,1,figsize=(w,h))
    
    inputDir = "results/{}/{}/".format(foldername,filename)
    
    # Load data
    try:
        inputs,outputs,profiles = readdata(inputDir)
    except:
        print('{} does not exist'.format(inputDir))
        return
    
    # Calculate bulk modulus from PREM
    bulk_modulus = premvp(profiles.r)**2*premdensity(profiles.r)
    
    # Calculate vp using slurry density and PREM bulk modulus
    vp_slurry = np.sqrt(bulk_modulus/profiles.density)
    
    # Calculate FVW P wave speed (Ohtaki et al. 2015, fig 11a)
    x = profiles.r/earth_radius
    vp_fvw = 3.3*x[0]-3.3*x +10.33
    
    # Look up AK135
    radius_ak135 = ak135radius()
    vp_ak135 = ak135vp(radius_ak135)
    
    # Check density
    # ax1.plot(profiles.r*1e-3,premdensity(profiles.r),'k--')
    # ax1.plot(profiles.r*1e-3,profiles.density)
    
    max_diff = np.max((vp_fvw-vp_slurry*1e-3)/vp_fvw*100)
    print('Maximum difference is {:.2f}%'.format(max_diff))
    
    # Plot P wave speed
    ax.plot(profiles.r*1e-3,vp_slurry*1e-3,color='darkgrey',lw=2,label='slurry') #(km/s)
    ax.vlines(profiles.r[0]*1e-3,vp_slurry[0]*1e-3,10.4,color='darkgrey',lw=2)
    ax.plot(profiles.r*1e-3,vp_fvw,color='blue',lw=2,ls=':',label='Ohtaki et al. (2015)')
    ax.vlines(profiles.r[0]*1e-3,vp_fvw[0],10.4,color='blue',lw=2,ls=':')
    ax.plot(profiles.r*1e-3,premvp(profiles.r)*1e-3,'k--',label='PREM')
    ax.vlines(profiles.r[0]*1e-3,premvp(profiles.r[0])*1e-3,10.4, 'k', linestyle='--')
    ax.plot(radius_ak135*1e-3,vp_ak135*1e-3,'k',label='ak135')
    ax.vlines(radius_ak135[0]*1e-3,vp_ak135[0]*1e-3,10.4, 'k')
    
    ax.legend(fontsize=11.5)
    ax.set(xlabel="Radius (km)")
    ax.set(ylabel="P wave speed (km/s)")
    ax.set_xlim([1200,profiles.r.iloc[-1]*1e-3])
    ax.set_ylim([10.18,10.4])
    plt.yticks(np.arange(10.2,10.4,0.1))
    
    if saveOn==1:
        saveDir='figures/seismic/'
        if not os.path.exists(saveDir+foldername):
            os.makedirs(saveDir+foldername)
        fig.savefig(saveDir+foldername+"/"+filename+".pdf",format='pdf', dpi=200, bbox_inches='tight')
        fig.savefig(saveDir+foldername+"/"+filename+".png",format='png', dpi=200, bbox_inches='tight')
        print('Figure saved as {}'.format(saveDir+foldername+"/"+filename+".pdf"))
    plt.show()