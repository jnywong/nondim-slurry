#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 17:31:41 2020

@author: jennywong

###############################################################################
# script_regime.py                                                            #
###############################################################################

# Plot regime diagram presented in Wong et al. (in prep)(see also Wong et al. 2018).

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
import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from slurpy.postprocess import regime

# %% MODEL INPUTS
layer_thickness = 250e3
thermal_conductivity = np.array([100.,30.])
h=0.5 # stepsize of heat flux through parameter space
# csb_heatfluxes=np.arange(5,6.05,h) # (TW)
# icb_heatfluxes=np.arange(2,4.05,h)
csb_heatfluxes=np.arange(0.05,15.05,h) # (TW)
icb_heatfluxes=np.arange(0.05,5.05,h) # (TW)

figAspect = 0.75
saveOn=0
xLower = 0
xUpper = 2500
yLower = 0
yUpper = 3
n_levels = 13

w, h = plt.figaspect(figAspect)
fig = plt.figure(constrained_layout=True, figsize=(1.5*w,1.5*h))
spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)
ax1 = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1], sharex = ax1, sharey = ax1)
ax3 = fig.add_subplot(spec[1, 0], sharex = ax1, sharey = ax1)
ax4 = fig.add_subplot(spec[1, 1], sharex = ax1, sharey = ax1)

fig2, (ax2a, ax2b) = plt.subplots(1,2,figsize=(3*w,1.5*h))
# Load data
for k in thermal_conductivity:
    if k ==100.:
        (Pe_100,St_100,density_jump_100,Q_cmb_100,stable_100) = \
                        regime(layer_thickness,k,csb_heatfluxes,icb_heatfluxes)
    elif k==30.:
        (Pe_30,St_30,density_jump_30,Q_cmb_30,stable_30) = \
                        regime(layer_thickness,k,csb_heatfluxes,icb_heatfluxes)

Z = Q_cmb_100*1e-12
levels = np.linspace(5,17,n_levels)
h1 = ax1.contourf(Pe_100,St_100,Z,levels,cmap='YlOrRd_r',extend='both')
c1 = fig.colorbar(h1, ax=ax2, shrink=1)
c1.ax.set_ylabel(r'CMB heat flux ($TW$)')

h2 = ax2.contourf(Pe_30,St_30,Q_cmb_30*1e-12,cmap='YlOrRd_r')

h3 = ax3.contourf(Pe_100,St_100,density_jump_100,cmap='YlGnBu_r')
c2 = fig.colorbar(h3, ax=ax4, shrink=1)
c2.ax.set_ylabel(r'density jump ($kgm^{-3}$)')

h4 = ax4.contourf(Pe_30,St_30,density_jump_30,cmap='YlGnBu_r')

# h5 = ax2a.contourf(Pe_100,St_100,stable_100,levels=[0,1,2,3]) #,color='k')
# ax2b.contourf(Pe_30,St_30,stable_30,levels=[0,1,2,3]) #,color='k')
# fig2.colorbar(h5, ax=ax2b, shrink=1)

ax2a.contourf(Pe_100,St_100,stable_100) #,levels=[0,1,2,3]) #,color='k')
ax2a.scatter(Pe_100,St_100,c=stable_100)

ax1.axhline(0.69)
ax1.set_xlim((xLower,xUpper))
ax1.set_ylim((yLower,yUpper))
ax2a.set_ylim((yLower,yUpper))
ax3.set_xlabel(r'$Pe$')
ax4.set_xlabel(r'$Pe$')
ax1.set_ylabel(r'$St$')    
ax3.set_ylabel(r'$St$')    

if saveOn==1:
        saveDir='figures/regime/'
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        fig.savefig(saveDir+"/regime.pdf",format='pdf', dpi=200, bbox_inches='tight')
        fig.savefig(saveDir+"/regime.png",format='png', dpi=200, bbox_inches='tight')
        print('Figure saved as {}'.format(saveDir+"/regime.pdf"))
plt.show()

