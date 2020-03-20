# -*- coding: utf-8 -*-
"""
Plot spherical profiles
"""

# Compare spherical profiles
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import cycler

from coreproperties import icb_radius, aO
from refparameters import getcsbmassoxygen,getL1,getL2,getPeclet,getStefan, \
    getLewis,getcsbradius,getfreezingspeed,getchangevolmelting
from readscripts import read_sphdata
from parametersearch import stablelayer
from lookup import premdensity

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

#------------------------------------------------------------------------------
# INPUTS
#------------------------------------------------------------------------------
saveOn=1
plotOn=0

#layer_thickness=np.array([150e3])
layer_thickness=np.array([150e3,200e3,250e3,300e3,350e3,400e3])
csb_heatflux=5 #3.5  
icb_heatflux=2.5 #1.5 
thermal_conductivity=100.0

mol_conc_oxygen_bulk=8. # light element concentration in the bulk of the OC (mol.%)
mol_conc_SSi = 8. # concentration of silicon/sulphur (mol.%)
self_diffusion= 0.98e-8

aspectRatio=0.75

# Axis limits
tempMin = 5400
tempMax = 5800
xiMin = 6
xiMax = 8
jMin = -3.5e-7
jMax = 0
denMin = 11900
denMax = 12200

#------------------------------------------------------------------------------

w,h= plt.figaspect(aspectRatio)*2
fig=plt.figure()
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,figsize=(w,h))
colors = plt.cm.tab10(np.linspace(0,1,layer_thickness.size)) 
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab10.colors)

if thermal_conductivity==100:
    folderName='highLe/'
else:
    folderName='lowLe/'

for i in (range(layer_thickness.size)):
    # %% Find dimensionless numbers
    csb_radius=getcsbradius(layer_thickness[i])
    mass_conc_O,acore=getcsbmassoxygen(mol_conc_oxygen_bulk,mol_conc_SSi)
    density0 = premdensity(csb_radius)
    deltaV_liquidFeO_solidFe=getchangevolmelting(mol_conc_oxygen_bulk,mol_conc_SSi,density0)
    freezing_speed=getfreezingspeed(icb_heatflux)
    
    L1,_,_=getL1(csb_radius)
    L2=getL2(mass_conc_O)
    St=getStefan(icb_heatflux,csb_heatflux,csb_radius)
    Le=getLewis(thermal_conductivity,self_diffusion,density0)
    Pe=getPeclet(freezing_speed,csb_radius,self_diffusion)
    
    # %% INPUT DIRECTORY
    str1=str(np.round(L1,2))
    str2=str(np.round(L2,2))
    str3=str(np.round(St,2))
    str4=str(np.round(Le,2))
    str5=str(np.round(Pe,2))
    
    fileName="L1" + str1 + "_L2" + str2 + "_Le" + str4 + "_St" + str3 + "_Pe" + str5 + "/" 
    inputDir="results/"+ folderName+ fileName
    
    data_in,data_out,data_profiles=read_sphdata(inputDir)
    
    density_jump = data_profiles.density.iloc[0]-data_profiles.density.iloc[-1]
    
    print('Density jump is {:.2f}'.format(density_jump))

    radius=(data_profiles['z'])*1e-3
    oxygen=data_profiles['oxygen']*acore/aO*100
    ax1.plot(radius,data_profiles['temp'],label='_nolegend_')#, color=colors[i])
    ax2.plot(radius,oxygen)#, color=colors[i])
    ax3.plot(radius,data_profiles['solidflux'])#,color=colors[i])
    ax4.plot(radius,data_profiles['density'])#,color=colors[i])
    
    stablelayer(icb_heatflux, \
                layer_thickness[i],csb_heatflux,thermal_conductivity,plotOn)

# Liquidus        
from lookup import liquidus
radius_liquidus=np.linspace(icb_radius,icb_radius+400e3)
temp_liquidus=liquidus(radius_liquidus)
ax1.plot(radius_liquidus*1e-3,temp_liquidus,'k--')

# Constant CSB oxygen
#ax2.hlines(y=8, xmin=radius[0], xmax=radius.iloc[-1],colors='k',linestyles='dashed')
       
# Zero solid flux
#ax3.hlines(y=0, xmin=radius[0], xmax=radius.iloc[-1],colors='k',linestyles='dashed')
 
# PREM
from lookup import premdensity
radius_prem=np.linspace(icb_radius,icb_radius+400e3)
density_prem=premdensity(radius_prem)
ax4.plot(radius_prem*1e-3,density_prem,'k--')
ax4.set(xlabel="Radius (km)",ylabel="Density ($\mathrm{kg m^{-3}}$)")

# Axis titles
ax1.set(ylabel="Temperature (K)")
ax2.set(ylabel="Oxygen (mol.%)")
ax3.set(xlabel="Radius (km)",ylabel="Solid flux ($\mathrm{kg m^{-2} s^{-1}}$)")


# Legend
ax1.legend(['Liquidus (Davies et al. 2015)'],fontsize=11,loc=1) #,
#           bbox_to_anchor=(1, 0.5))
ax4.legend(['150 km','200 km','250 km','300 km','350 km','400 km','PREM'],
           fontsize=11,loc=1) #,
#           bbox_to_anchor=(1, 0.5))
#ax4.legend([r'233 $kg m^{-3}$',r'280 $kg m^{-3}$',r'311 $kg m^{-3}$',r'328 $kg m^{-3}$',r'338 $kg m^{-3}$',r'342 $kg m^{-3}$','PREM'],fontsize=10,loc='center left',
#           bbox_to_anchor=(1, 0.5))

# Axis limits
ax1.set_ylim([tempMin,tempMax])
ax2.set_ylim([xiMin,xiMax])
ax3.set_ylim([jMin,jMax])
ax4.set_ylim([denMin,denMax])

# Compare to poster
#ax1.set_ylim([5400,5750])
#ax2.set_ylim([6.8,8])
#ax3.set_ylim([-5e-7,0])
#ax4.set_ylim([11950,12200])

ax1.set_xlim([1220,(icb_radius+400e3)*1e-3])

# Subfigure labels
#ax1.text(1225,tempMax-23,'(a)',fontsize=14)
#ax2.text(1225,xiMax - 0.12,'(b)',fontsize=14)
#ax3.text(1225,jMax - .2e-7,'(c)',fontsize=14)
#ax4.text(1225,denMax - 20,'(d)',fontsize=14)
 
# Subfigure titles
ax1.set_title('(a)',x=0.95,y=1,fontsize=14)
ax2.set_title('(b)',x=0.95,y=1,fontsize=14)
ax3.set_title('(c)',x=0.95,y=1,fontsize=14)
ax4.set_title('(d)',x=0.95,y=1,fontsize=14)

plt.tight_layout()

if saveOn==1:
    if not os.path.exists('figures/profiles'):
        os.makedirs('figures/profiles')
#    saveName='icbhf'+str(np.round(icb_heatflux,2))+'_csbhf'+str(np.round(csb_heatflux,2)) + '_k'+str(np.round(thermal_conductivity,2)).replace(".","_")
    saveName='icbhf{}_csbhf{}_k{}'.format(str(np.round(icb_heatflux,2)).replace(".","_"),
                   str(np.round(csb_heatflux,2)).replace(".","_"),
                   str(np.round(thermal_conductivity,2)).replace(".","_"))
    plt.savefig('figures/profiles/'+saveName+'.pdf',format='pdf',dpi=200)
    
plt.show()