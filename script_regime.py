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
import matplotlib.cm as cm

from slurpy.postprocess import regime

# %% MODEL INPUTS
layer_thickness = 250e3
thermal_conductivity = np.array([100.,30.])
h=0.25 # stepsize of heat flux through parameter space
# csb_heatfluxes=np.arange(5,6.05,h) # (TW)
# icb_heatfluxes=np.arange(2,4.05,h)
csb_heatfluxes=np.arange(0.05,15.05,h) # (TW)
icb_heatfluxes=np.arange(0.05,5.05,h) # (TW)

saveOn=1
figAspect = 0.75
xLower = 0
xUpper = 2500
yLower = 0
yUpper = 3
n_levels = 13

dmin = 50.
dmax = 350.

w, h = plt.figaspect(figAspect)
fig = plt.figure(constrained_layout=True, figsize=(1.5*w,1.5*h))
spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)
nFig = 4
ax = []
ax.append(fig.add_subplot(spec[0, 0]))
ax.append(fig.add_subplot(spec[0, 1], sharex = ax[0], sharey = ax[0]))
ax.append(fig.add_subplot(spec[1, 0], sharex = ax[0], sharey = ax[0]))
ax.append(fig.add_subplot(spec[1, 1], sharex = ax[0], sharey = ax[0]))

fig2, (ax2a, ax2b) = plt.subplots(1,2,figsize=(3*w,1.5*h))
# Load data
for k in thermal_conductivity:
    if k ==100.:
        (Pe_100,St_100,density_jump_100,Q_cmb_100,stable_100) = \
                regime(layer_thickness,k,csb_heatfluxes,icb_heatfluxes,yUpper)
    elif k==30.:
        (Pe_30,St_30,density_jump_30,Q_cmb_30,stable_30) = \
                regime(layer_thickness,k,csb_heatfluxes,icb_heatfluxes,yUpper)

Z = Q_cmb_100*1e-12
# Z = ma.masked_where(Pe_100>xUpper+500,Z)
# Z = ma.masked_where(St_100>yUpper+2.5,Z)
levels = np.linspace(5,17,n_levels)
h1 = ax[0].contourf(Pe_100,St_100,Z,levels,cmap='YlOrRd_r',extend='both')
c1 = fig.colorbar(h1, ax=ax[1], shrink=1)
c1.ax.set_ylabel(r'CMB heat flux ($TW$)')

Z = Q_cmb_30*1e-12
# Z = ma.masked_where(Pe_30>xUpper+500,Z)
# Z = ma.masked_where(St_30>yUpper+2.5,Z)
# Z = ma.masked_where(Z>17,Z)
# Z = ma.masked_where(Z<5,Z)
h2 = ax[1].contourf(Pe_30,St_30,Z,levels,cmap='YlOrRd_r',extend='both') #,
                  # vmin = 5, vmax = 17)

Z = density_jump_100
levels = np.linspace(dmin,dmax,n_levels)
h3 = ax[2].contourf(Pe_100,St_100,Z,levels,cmap=cm.YlGnBu_r,extend = 'both')
# ax[2].scatter(Pe_100,St_100,c=Z)

c2 = fig.colorbar(h3, ax=ax[3], shrink=1)
c2.ax.set_ylabel(r'density jump ($kgm^{-3}$)')

h4 = ax[3].contourf(Pe_30,St_30,density_jump_30,levels,cmap='YlGnBu_r',extend = 'both')

# h5 = ax2a.contourf(Pe_100,St_100,stable_100,levels=[0,1,2,3]) #,color='k')
# ax2b.contourf(Pe_30,St_30,stable_30,levels=[0,1,2,3]) #,color='k')
# fig2.colorbar(h5, ax=ax2b, shrink=1)

# Stable

cOut = ax[0].contour(Pe_100,St_100,stable_100,levels=[0,1,2,3], alpha = 0.5, cmap = 'gray') #colors='k')
ax[1].contour(Pe_30,St_30,stable_30,levels=[0,1,2,3], alpha = 0.5, cmap = 'gray') #colors='k')
ax2a.contour(Pe_100,St_100,stable_100,levels=[0,1,2,3]) #,colors='k')
ax2a.scatter(Pe_100,St_100,c=stable_100)
ax2b.contour(Pe_30,St_30,stable_30,levels=[0,1,2,3]) #,colors='k')
# ax2b.contour(Pe_30,St_30,stable_30) #,levels=[0,1,2,3]) #,color='k')
ax2b.scatter(Pe_30,St_30,c=stable_30)
ax2a.set_xlim((xLower,xUpper))
ax2a.set_ylim((yLower,yUpper))
ax2b.set_xlim((xLower,xUpper))
ax2b.set_ylim((yLower,yUpper))

# Finishing touches
for axes in ax:
    axes.set_xlim((xLower,xUpper))
    axes.set_ylim((yLower,yUpper))
    axes.axhline(0.69,color='k')

ax[2].set_xlabel(r'$Pe$')
ax[3].set_xlabel(r'$Pe$')
ax[0].set_ylabel(r'$St$')    
ax[2].set_ylabel(r'$St$')    

if saveOn==1:
    saveDir='figures/regime/'
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    fig.savefig(saveDir+"regime.pdf",format='pdf', dpi=200, bbox_inches='tight')
    fig.savefig(saveDir+"regime.png",format='png', dpi=200, bbox_inches='tight')
    print('Figure saved as {}'.format(saveDir+"regime.pdf"))
plt.show()

