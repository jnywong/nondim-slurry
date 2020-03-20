#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 13:30:27 2019

@author: jennywong
"""

import os
import numpy as np
import matplotlib
#matplotlib.use('PS')
import matplotlib.pyplot as plt

from parametersearch import stablelayer, getphaseboundary
from savescripts import savephaseboundaries

#------------------------------------------------------------------------------
# REGIMEDIAGRAM_ALL.PY
#------------------------------------------------------------------------------
# Reads all the data to categorise solutions with:
#   - stable density
#   - unstable density
#   - unstable oxygen profile but stable density
#   - unstbale oxygen and unstable density

# PARAMETERSEARCH ASSUMES D=150km FOR PARTIAL STABILITY

#------------------------------------------------------------------------------
# SETTINGS
#------------------------------------------------------------------------------

saveOn=1
plotOn=0

# Select parameters
layer_thickness=np.array([150e3]) #(km)
#csb_heatfluxes=np.array([5]) # (TW)
thermal_conductivity=np.array([30.0]) # (W m^-1 K^-1)
#icb_heatfluxes=np.array([2.5]) # (TW)

h=0.05
csb_heatfluxes=np.arange(0.1,20.05,h) # (TW) 
icb_heatfluxes=np.arange(0.1,5.55,h) # (TW)
switch=2 # 0=box single eig., 1=box two eig., 2=spherical two eig.

#------------------------------------------------------------------------------

# Pre-allocate arrays
nPe=icb_heatfluxes.size
nSt=csb_heatfluxes.size
St=np.zeros([nPe,nSt])
Pe=np.zeros([nPe,nSt])
density_jump=np.zeros([nPe,nSt])
cmb_heatflux=np.zeros([nPe,nSt])
stable_density=np.zeros([nPe,nSt])
stable_oxygen=np.zeros([nPe,nSt])
solution=np.zeros([nPe,nSt])
density_exceeds_prem=np.zeros([nPe,nSt])
colors=np.zeros([nPe,nSt])
markers=np.zeros([nPe,nSt])
delta=np.zeros([nPe,nSt])
F=np.zeros([nPe,nSt])
icb_speed=np.zeros([nPe,nSt])

# PLOT
width, height = plt.figaspect(0.75)
fig1=plt.figure
fig1,ax1=plt.subplots(1,1,figsize=(width,height))

Pe_stable=[]
St_stable=[]
Pe_unstable=[]
St_unstable=[]
Pe_lower=[]
St_lower=[]
Pe_upper=[]
St_upper=[]

for i in range(nPe):
    for j in range(nSt):
        density_jump[i,j],cmb_heatflux[i,j],stable_density[i,j],stable_oxygen[i,j], \
                solution[i,j],density_exceeds_prem[i,j],L1,L2,St[i,j],Pe[i,j],Le, \
                F[i,j], icb_speed[i,j] =stablelayer(icb_heatfluxes[i], \
                layer_thickness,csb_heatfluxes[j],thermal_conductivity,plotOn)
        if solution[i,j]==1:
            if density_exceeds_prem[i,j]==0:
                delta[i,j]=0
            elif density_exceeds_prem[i,j]==1:            
                # Scatter plot
                if (stable_density[i,j]==1) & (stable_oxygen[i,j]==1):
                    ax1.scatter(Pe[i,j],St[i,j],c='#bfedff',marker='s') # stable density and stable O blue
                    delta[i,j]=1
                elif (stable_density[i,j]==1) & (stable_oxygen[i,j]==0):
                    ax1.scatter(Pe[i,j],St[i,j],c='#f0f06c',marker='s') # stable density unstable O yellow
                    delta[i,j]=2
                elif (stable_density[i,j]==2) & (stable_oxygen[i,j]== 1):
                    ax1.scatter(Pe[i,j],St[i,j],c='#f0c665',marker='s') # partially stratified stable O orange
                    delta[i,j]=3                 
                elif (stable_density[i,j]==2) & (stable_oxygen[i,j]== 0):
                    ax1.scatter(Pe[i,j],St[i,j],c='#f0c665',marker='s') # partially stratified unstable O orange
                    delta[i,j]=4 
                elif (stable_density[i,j]==0):
                    ax1.scatter(Pe[i,j],St[i,j],c='#d46666',marker='s') # unstable density red
                    delta[i,j]=0
                # Line plot
                if (stable_density[i,j]==1) and (stable_oxygen[i,j]==0) and (stable_oxygen[i-1,j]==1) and (i>0):
                    Pe_stable.append(Pe[i,j])
                    St_stable.append(St[i,j])
                if (stable_density[i,j]==1) and (stable_density[i-1,j]==0) and (i>0) and (St[i,j]<4):
                    Pe_lower.append(Pe[i,j])
                    St_lower.append(St[i,j])    
                if (stable_density[i,j]==1) and (density_exceeds_prem[i-1,j]==0) and (solution[i-1,j]==1) and (i>0): #and (St<6)
                    Pe_unstable.append(Pe[i,j])
                    St_unstable.append(St[i,j])    
        elif solution[i,j]==0:
            ax1.scatter(Pe[i,j],St[i,j],c='black',marker='s')
            delta[i,j]=5
            if (stable_density[i-1,j]==1) and (i>0):
                Pe_upper.append(Pe[i,j])
                St_upper.append(St[i,j]) 
                
ax1.set(yscale='log')

if saveOn==1:
    if not os.path.exists('figures/regime_diagram'):
        os.makedirs('figures/regime_diagram')
    saveName='sph_L1'+str(np.round(L1,2))+'_L2'+str(np.round(L2,2)) + '_Le'+str(np.round(Le,2))
    plt.savefig('figures/regime_diagram/'+saveName+'_raw_boundary.png',format='png',dpi=200)
    
    if not os.path.exists('phaseboundaries/'+saveName+'/'):
        os.makedirs('phaseboundaries/'+saveName+'/')
    outputDir='phaseboundaries/'+saveName+'/'
    savephaseboundaries(outputDir,nPe,nSt,Pe,St,density_jump,cmb_heatflux,delta,
                        F,icb_speed)
    
#plt.show()
