#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 15:22:04 2020

@author: wong
"""
import pickle
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os

from matplotlib.ticker import PercentFormatter

from lookup import premdensity

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

#------------------------------------------------------------------------------
# INPUTS
#------------------------------------------------------------------------------
csb_oxygen = np.arange(2,12.5,0.5)
csb_oxygen0 = 8.0

csb_temp = np.arange(4500,6100,100)
csb_temp0 = 5547

aspectRatio=0.75
max_den = 12250

saveOn=1
#------------------------------------------------------------------------------

#%% Oxygen
w, h = plt.figaspect(aspectRatio)
fig1, (ax1,ax2) = plt.subplots(1,2,figsize=(2*w,h),sharey=True)
nOxy = csb_oxygen.size
colors=plt.cm.GnBu(np.linspace(0.4,1,nOxy))
den_jump=[]

for i in range(nOxy):
    filename = 'sensitivity/xi_{:.1f}'.format(csb_oxygen[i]).replace('.','_')
    with open(filename, 'rb') as f:
        (radius,temp,xi,solidFlux,density)=pickle.load(f)
    if i ==0 or i == nOxy-1:
        ax1.plot(radius*1e-3,density,color=colors[i], linewidth = 2, 
                 label =r'$\xi_O=${:.1f} mol.%'.format(csb_oxygen[i]))
    # Reference case        
    elif csb_oxygen[i]==csb_oxygen0:
        den_jump0 = density[0]-density[-1]
        ax1.plot(radius*1e-3,density,color='silver', linewidth = 2,
                 label=r'$\xi_O=$8.0 mol.%')
    else:        
        ax1.plot(radius*1e-3,density,color=colors[i])    
    den_jump.append(density[0]-density[-1])
   
# PREM
density_prem=premdensity(radius) 
ax1.plot(radius*1e-3,density_prem, 'k', linestyle = '--')

ax1.legend(fontsize=11.5)
ax1.set_xlim([radius[0]*1e-3,radius[-1]*1e-3])
ax1.set_ylim([12080,max_den])
ax1.set(xlabel="Radius (km)",ylabel="Density ($\mathrm{kg m^{-3}}$)")
#plt.xticks(np.arange(1240,1370,20))

den_low = (den_jump[0]-den_jump0)/den_jump0*100
den_high = (den_jump[-1]-den_jump0)/den_jump0*100
print('Oxygen: Density jump ranges from {:.2f}% to {:.2f}% of reference'.format(den_low, den_high))

# Oxygen versus density jump
fig3, ax3 = plt.subplots(1,1,figsize=(w,h))
ax3.yaxis.set_major_formatter(PercentFormatter())
y = (den_jump - den_jump0)/den_jump0*100
ax3.plot(csb_oxygen,y,color=colors[-1],lw=3)
#idx = int(np.argwhere(csb_oxygen == csb_oxygen0))
#ax3.plot(csb_oxygen0,den_jump[idx],'*')
#ax3.set_xlim([0,12.5])
#ax3.set(xlabel="CSB oxygen (mol.%)",ylabel="Change in density jump")

#%% Temperature
#fig2, ax2 = plt.subplots(1,1,figsize=(w,h))
csb_temp = np.append(csb_temp,csb_temp0)
csb_temp = np.sort(csb_temp)
nTemp = csb_temp.size
colors=plt.cm.OrRd(np.linspace(0.4,1,nTemp))
den_jump=[]

for i in range(nTemp):
    filename = 'sensitivity/temp_{:.0f}'.format(csb_temp[i]).replace('.','_')
    with open(filename, 'rb') as f:
        (radius,temp,xi,solidFlux,density)=pickle.load(f)
    if i ==0 or i == nTemp-1:
        ax2.plot(radius*1e-3,density,color=colors[i], linewidth = 2, 
                 label =r'$T_l=${:.0f} K'.format(csb_temp[i]))
    # Reference case        
    elif csb_temp[i]==csb_temp0:
        den_jump0 = density[0]-density[-1]
        ax2.plot(radius*1e-3,density,color='silver', linewidth = 2,
                 label=r'$T_l=$5457 K')
    else:        
        ax2.plot(radius*1e-3,density,color=colors[i])    
    den_jump.append(density[0]-density[-1])
   
# PREM
density_prem=premdensity(radius) 
ax2.plot(radius*1e-3,density_prem, 'k', linestyle = '--', label = r'PREM')

ax2.legend(fontsize=11.5)
ax2.set_xlim([radius[0]*1e-3,radius[-1]*1e-3])
ax2.set(xlabel="Radius (km)") #,ylabel="Density ($\mathrm{kg m^{-3}}$)")
#ax2.set_xticks(np.arange(1225,1380,25))

den_low = (den_jump[0]-den_jump0)/den_jump0*100
den_high = (den_jump[-1]-den_jump0)/den_jump0*100
print('Temperature: Density jump ranges from {:.2f}% to {:.2f}% of reference'.format(den_low, den_high))

# Temperature vs density jump
#ax4 = ax3.twiny()
#y = (den_jump - den_jump0)/den_jump0*100
#ax4.plot(csb_temp,y,color=colors[-1],lw=3)
#idx = int(np.argwhere(csb_temp == csb_temp0))
#ax4.plot(csb_temp0,y[idx],'*', color = 'yellow', ms = 18, markeredgecolor = 'k')
#ax4.set_xlim([4400,6200])
#ax4.set_ylim([-35,35])
#ax4.set(xlabel="CSB temperature (K)",ylabel="Change in density jump")

if saveOn==1:
    saveDir='figures/sensitivity/'
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    fig1.savefig(saveDir+"sensitivity.pdf",format='pdf', dpi=200, bbox_inches='tight')
#    fig2.savefig(saveDir+"temp.pdf",format='pdf', dpi=200, bbox_inches='tight')  
    fig3.savefig(saveDir+"compare.pdf",format='pdf', dpi=200, bbox_inches='tight') 