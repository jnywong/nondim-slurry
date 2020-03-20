#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 13:53:38 2019

@author: wong
"""

#------------------------------------------------------------------------------
# SEDIMENTATION.PY
# Sensitivity of sedimentation coefficient to free parameter N
#------------------------------------------------------------------------------
import numpy as np
import socket
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import os

from matplotlib.pyplot import cm

from parameter_range import getdimensionless
from readscripts import read_sphdata
from postprocess import slurrydensity
from lookup import premdensity
from coreproperties import icb_radius

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

#%% PLEASE ENTER INPUTS:
# Input parameters
layer_thickness = 150e3
icb_heatflux=2.5
csb_heatflux=5
thermal_conductivity=100.

# Range of sedimentation constants
sed_lower=-5
sed_upper=-1

# Figure aspect ratio
figAspect=0.75

# Save figures?
saveOn=1
#%% Load data
(L1,L2,Pe,St,Le)=getdimensionless(layer_thickness,icb_heatflux,
                                        csb_heatflux,thermal_conductivity)
if Le>400:
    subDir='highLe'
else:
    subDir='lowLe'
str1=str(np.round(L1,2))
str2=str(np.round(L2,2))
str3=str(np.round(Le,2))
str4=str(np.round(St,2))
str5=str(np.round(Pe,2))
filename="{}/L1{}_L2{}_Le{}_St{}_Pe{}/".format(subDir,str1,str2,str3,str4,str5)

if socket.gethostname()=='JenBook.local':
    directory="./"
elif socket.gethostname()=='imacjenny.ipgp.fr':
    directory="./"
    
(input_data,output_data,profile_data)=read_sphdata(directory+"results/"+filename)

# Convert from pandas to numpy arrays
radius=profile_data.z.to_numpy()
temp=profile_data.temp.to_numpy()
xi=profile_data.oxygen.to_numpy()
solidflux=profile_data.solidflux.to_numpy()
mol_conc_oxygen_bulk=input_data.oxygen_bulk.to_numpy()
mol_conc_SSi=input_data.siliconSulphur_bulk.to_numpy()
#sedimentation_constant=input_data.sedimentationConstant.to_numpy()

#%% Calculate slurry density
nSed=sed_upper-sed_lower+1
sedimentation_constant=np.logspace(sed_lower,sed_upper,num=nSed)
colors=plt.cm.copper_r(np.linspace(0.4,1,nSed))

w, h = plt.figaspect(figAspect)
gs = gridspec.GridSpec(1, 1, height_ratios=[1])
fig, (ax1,ax2) = plt.subplots(2,1,sharex=True,figsize=(w,2*h))

# Format legend into scientific notation
f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
g = lambda x,pos : "{}".format(f._formatSciNotation('%1.10e' % x))
fmt = mticker.FuncFormatter(g)

flucMax=np.zeros((nSed))
flucMin=np.zeros((nSed))

for i in range(nSed):
    (density,phi,temp_fluc,xi_fluc,phi_fluc,density_fluc)=slurrydensity(radius,temp,xi,
        solidflux,layer_thickness,mol_conc_oxygen_bulk,mol_conc_SSi,
        sedimentation_constant[i])
    flucMax[i]=max(density_fluc.max(),xi_fluc.max())
    flucMin[i]=temp_fluc.min()
#    c=next(color)
    ax1.plot(radius*1e-3,density,label=r"$k_\phi={}$".format(fmt(sedimentation_constant[i])),
             color=colors[i])
    ax2.plot(radius*1e-3,phi,color=colors[i]) #label=r"$k_\phi={}$".format(fmt(sedimentation_constant[i])))
             
    
density_prem=premdensity(radius)
ax1.plot(radius*1e-3,density_prem,'k--', label='PREM')
ax1.set(ylabel="Density ($\mathrm{kg m^{-3}}$)") #,yscale="log")
ax2.set(xlabel="Radius (km)",ylabel="Solid fraction",yscale='log')
#ax2.axhline(0.6,color='k',linestyle='--') # rheological transition
fig.legend(loc='center right', bbox_to_anchor=(1.4, 0.5),fontsize = 11.5)
ax1.set_xlim([radius[0]*1e-3,radius[-1]*1e-3])
ax2.set_ylim([1e-4,1])
#plt.yscale('log')

ax1.set_title('(a)',x=0.95,y=1,fontsize=14)
ax2.set_title('(b)',x=0.95,y=1,fontsize=14)

if saveOn==1:
    saveDir=directory+'figures/sedimentation/'
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    fig.savefig(saveDir+"sedimentation.pdf",format='pdf', dpi=200, bbox_inches='tight')

plt.show()
