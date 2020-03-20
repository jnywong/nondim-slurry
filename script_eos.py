#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 17:29:09 2020

@author: wong
"""

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from scipy.integrate import cumtrapz

from readscripts import read_sphdata
from parameter_range import getdimensionless
from refparameters import getfreezingspeed, getcsbmassoxygen, getchangevolmelting
from lookup import premdensity, premgravity, liquidus
from coreproperties import latent_heat, heat_capacity, icb_radius, \
    deltaV_solidFe_liquidFe, density_solidFe, aO
from mainroutines import slurrydensity


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)
#------------------------------------------------------------------------------
aspectRatio=0.75

#saveNames = ['L10.16_L20.02_Le1180.89_St1.59_Pe1959.06','L10.22_L20.02_Le1195.97_St1.13_Pe2316.16']
saveNames = ['L10.16_L20.02_Le1180.89_St0.95_Pe1959.06','L10.16_L20.02_Le1180.89_St1.59_Pe1959.06']

w, h = plt.figaspect(aspectRatio)
fig1, ax1a = plt.subplots(1,1,figsize=(w,h))
fig2, (ax2a,ax2b) = plt.subplots(1,2,figsize=(2*w,h))
fig3, (ax3a,ax3b) = plt.subplots(1,2,figsize=(2*w,h))
fig4, (ax4a,ax4b) = plt.subplots(1,2,figsize=(2*w,h))
fig5, (ax5a,ax5b) = plt.subplots(1,2,figsize=(2*w,h))

ax1a.axhline(0,color='k', linestyle='--')
colors=plt.cm.tab10(np.linspace(0,1,6))
for i in range(len(saveNames)):
    saveName = saveNames[i]
    outputDir = 'results/highLe/{}/'.format(saveName)
    inputs,outputs,profiles = read_sphdata(outputDir)
    
    layer_thickness = pd.to_numeric(inputs.layerThickness.iloc[0])
    icb_heatflux = pd.to_numeric(inputs.icb_heatflux.iloc[0])
    csb_heatflux = pd.to_numeric(inputs.csb_heatflux.iloc[0])
    thermal_conductivity = pd.to_numeric(inputs.thermalConductivity.iloc[0])
    
    csb_radius = icb_radius + layer_thickness
    freezing_speed=getfreezingspeed(icb_heatflux)
    rho0 = premdensity(csb_radius)
    scale_temp = csb_heatflux*1e12/(4*np.pi*csb_radius**2*rho0\
                                *heat_capacity*freezing_speed)
    scale_xi,acore = getcsbmassoxygen(pd.to_numeric(inputs.oxygen_bulk)[0],
                                      pd.to_numeric(inputs.siliconSulphur_bulk)[0])
    scale_j = density_solidFe*freezing_speed
    
    
    (L1,L2,Pe,St,Le) = getdimensionless(layer_thickness,icb_heatflux,csb_heatflux,
                         thermal_conductivity)

    print('{}: F = {:.2f}'.format(saveName,pd.to_numeric(outputs.F[0])))
    v_out =(outputs.snowSpeed + freezing_speed)/freezing_speed
    print(r'{}: v  = {:.2f}'.format(saveName,pd.to_numeric(v_out[0])))    

    # Plot
    x = pd.to_numeric(profiles.z)/csb_radius
    temp = pd.to_numeric(profiles.temp)/scale_temp
    xi = pd.to_numeric(profiles.oxygen)/scale_xi
    solid_flux = pd.to_numeric(profiles.solidflux)/scale_j
    dtemp = pd.to_numeric(profiles.temp_grad)*csb_radius/scale_temp
    dxi = pd.to_numeric(profiles.xi_grad)*csb_radius/scale_xi
    dj = pd.to_numeric(profiles.j_grad)*csb_radius/scale_j 
    # Plot density
    density,_,temp_fluc,xi_fluc,phi_fluc,density_fluc = slurrydensity(profiles.z.to_numpy(),
                                                 profiles.temp.to_numpy(),
                                                 profiles.oxygen.to_numpy(),
                                                 profiles.solidflux.to_numpy(),
                                                 layer_thickness,
                                                 8, 8, 1e-2) # numeric

    drho_fluc = np.gradient(density_fluc,x)
    dtemp_fluc = np.gradient(temp_fluc,x)
    dxi_fluc = np.gradient(xi_fluc,x)
    dphi_fluc = np.gradient(phi_fluc,x)
    if i==0:
        ax1a.plot(x*csb_radius*1e-3,drho_fluc,color='darkgrey',lw=4,ls=':')        
        ax1a.plot(x*csb_radius*1e-3,dtemp_fluc,color='darkred',lw=2,ls=':')
        ax1a.plot(x*csb_radius*1e-3,dxi_fluc,color='royalblue',lw=2,ls=':')
        ax1a.plot(x*csb_radius*1e-3,dphi_fluc,color='peru',lw=2,ls=':')
#        ax1a.axvline(csb_radius*1e-3,color=colors[i])
        ax1a.axvspan(x.iloc[0]*csb_radius*1e-3, csb_radius*1e-3, facecolor='tab:blue', alpha=0.1)
        ax2a.plot(x*csb_radius*1e-3,density_fluc)
        ax2b.plot(x*csb_radius*1e-3,density-density_fluc,ls=':',label = 'hydrostatic')
        ax2b.plot(x*csb_radius*1e-3,density,ls=':',label = 'total')
        ax3a.plot(x,temp)
        ax3b.plot(x,dtemp)
        ax4a.plot(x,xi)
        ax4b.plot(x,dxi)
        ax5a.plot(x,solid_flux)
        ax5b.plot(x,dj)
    else:
        ax1a.plot(x*csb_radius*1e-3,drho_fluc,color='darkgrey',lw=4,label=r'$\mathrm{d} \rho^\prime = \rho_{sl} (-\alpha_T \mathrm{d}T^\prime -  \alpha_\xi \mathrm{d}\xi^\prime +  (\alpha_\phi + \alpha_\xi \xi)\mathrm{d}\phi^\prime)$')
        ax1a.plot(x*csb_radius*1e-3,dtemp_fluc,color='darkred',lw=2,label=r'$-\rho_{sl} \alpha_T \mathrm{d}T^\prime$')
        ax1a.plot(x*csb_radius*1e-3,dxi_fluc,color='royalblue',lw=2,label=r'$-\rho_{sl} \alpha_\xi \mathrm{d}\xi^\prime$')
        ax1a.plot(x*csb_radius*1e-3,dphi_fluc,color='peru',lw=2,label=r'$\rho_{sl} (\alpha_\phi + \alpha_\xi \xi)\mathrm{d}\phi^\prime$')
        ax1a.axvspan((icb_radius+150e3)*1e-3, csb_radius*1e-3, facecolor='tab:brown', alpha=0.1)
        ax2a.plot(x*csb_radius*1e-3,density_fluc)
        ax2b.plot(x*csb_radius*1e-3,density-density_fluc,label = 'hydrostatic')
        ax2b.plot(x*csb_radius*1e-3,density, label = 'total')
        ax3a.plot(x,temp)
        ax3b.plot(x,dtemp)
        ax4a.plot(x,xi)
        ax4b.plot(x,dxi)
        ax5a.plot(x,solid_flux)
        ax5b.plot(x,dj)
#    axs[i].xlabel('Radius (km)')

ax1a.legend(fontsize=7.6,loc=4) #, bbox_to_anchor=(0.5,-0.45))
ax1a.set(xlabel=r'Radius (km)')
#ax1a.set(ylabel=r'$\frac{\mathrm{d} \rho^\prime}{\mathrm{d}r} \ (\mathrm{kg m^{-2}}$)')
ax1a.set_xlim([x.iloc[0]*csb_radius*1e-3,x.iloc[-1]*csb_radius*1e-3])
ax1a.set_ylim([-1500,500])

ax1b = ax1a.twiny()
ax1b.set_xticks(np.array([(icb_radius-612e3)*1e-3,csb_radius*1e-3]))
labels=['150','400']
ax1b.set_xticklabels(labels)
ax1b.set_xlabel('Layer thickness (km)')


ax2a.set(xlabel=r'Radius (km)',ylabel = r'$\rho^\prime$')
ax2b.set(xlabel=r'Radius (km)',ylabel= r'$\rho$')
ax2b.legend()
fig2.tight_layout()

ax3a.set(xlabel=r'Radius (km)',ylabel = r'$T$')
ax3b.set(xlabel=r'Radius (km)',ylabel= r'$\frac{dT}{dr}$')
#ax3b.legend()
fig3.tight_layout()

ax4a.set(xlabel=r'Radius (km)',ylabel = r'$\xi$')
ax4b.set(xlabel=r'Radius (km)',ylabel= r'$\frac{d \xi}{dr}$')
#ax4b.legend()
fig4.tight_layout()

ax5a.set(xlabel=r'Radius (km)',ylabel = r'$j$')
ax5b.set(xlabel=r'Radius (km)',ylabel= r'$\frac{d j}{dr}$')
#ax5b.legend()
fig5.tight_layout()

fig1.savefig('figures/eos/drho.pdf',format='pdf',dpi=200, bbox_inches='tight') 