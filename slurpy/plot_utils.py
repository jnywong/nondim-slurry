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
from matplotlib import cm
import pickle
import os
import matplotlib as mpl
mpl.rcParams['text.usetex'] = False
import matplotlib.gridspec as gridspec

from slurpy.data_utils import readdata, get_outputDir
from slurpy.getparameters import getcsbmassoxygen, getKphi, getphi
from slurpy.coreproperties import icb_radius, earth_radius, aO, density_solidFe
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
    ax4.set(xlabel="Radius (km)",ylabel="Density ($\mathrm{kg m^{-3}}$)")

    # Liquidus
    radius_liquidus=np.linspace(icb_radius,csb_radius)
    temp_liquidus=liquidus(radius_liquidus)
    ax1.plot(radius_liquidus*1e-3,temp_liquidus,'k--', label = 'Liquidus (Davies et al. 2015)')
    ax1.legend()

    # Seismology
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
                     label =r'$T^{{sl}}=${:.0f} K'.format(csb_temp[i]))
        # Reference case
        elif csb_temp[i]==csb_temp0:
            den_jump0 = density[0]-density[-1]
            ax1.plot(radius*1e-3,density,color='silver', linewidth = 2,
                     label=r'$T^\mathregular{{sl}}=$5457 K')
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
                     label =r'$\xi^{{sl}}=${:.1f} mol.%'.format(csb_oxygen[i]))
        # Reference case
        elif csb_oxygen[i]==csb_oxy0:
            den_jump0 = density[0]-density[-1]
            ax2.plot(radius*1e-3,density,color='silver', linewidth = 2,
                     label=r'$\xi^{{sl}}=$8.0 mol.%')
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
    colors=plt.cm.gnuplot_r(np.linspace(0.4,1,nSed))

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
def plot_seismic(layer_thickness, thermal_conductivity,
                             icb_heatflux, csb_heatflux,saveOn,figAspect=0.75):
    w, h = plt.figaspect(figAspect)
    fig, ax = plt.subplots(1,1,figsize=(1.25*w,1.25*h))

    n = layer_thickness.size*thermal_conductivity.size \
        *icb_heatflux.size*csb_heatflux.size

    if n!=1:
        fig_label = [r'high $Le$', r'low $Le$']
        start = 0.2
        stop = 0.5
        cm_subsection = np.linspace(start, stop, n)
        colors = [ cm.gray_r(x) for x in cm_subsection ]

    # Load data
    k=0
    for w,x,y,z in [(w,x,y,z) for w in layer_thickness for x in icb_heatflux for y in csb_heatflux for z in thermal_conductivity]:
        foldername, filename = get_outputDir(w,x,y,z)
        inputDir = "results/{}/{}/".format(foldername,filename)

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

        max_diff = np.max((vp_fvw-vp_slurry*1e-3)/vp_fvw*100)
        print('Maximum difference with Ohtaki et al. (2015) is {:.2f}%'.format(max_diff))
        max_diff = np.max((premvp(profiles.r)-vp_slurry)/premvp(profiles.r)*100)
        print('Maximum difference with PREM is {:.2f}%'.format(max_diff))
        print('Density on slurry side of ICB is {:.2f}'.format(profiles.density[0]))
        density_jump = profiles.density[0] - profiles.density.iloc[-1]
        print('Density jump is {:.2f}'.format(density_jump))
        rho_bod = density_solidFe - profiles.density[0]
        print('Delta rho bod is {:.2f}'.format(rho_bod))
        rho_mod = rho_bod + density_jump
        print('Delta rho mod is {:.2f}'.format(rho_mod))

        # Plot P wave speed
        if n==1:
            ax.plot(profiles.r*1e-3,vp_slurry*1e-3,color='darkgrey',lw=2,label='slurry') #(km/s)
            ax.vlines(profiles.r[0]*1e-3,vp_slurry[0]*1e-3,10.4,color='darkgrey',lw=2)
        else:
            ax.plot(profiles.r*1e-3,vp_slurry*1e-3,color=colors[k],lw=2,label=fig_label[k]) #(km/s)
            ax.vlines(profiles.r[0]*1e-3,vp_slurry[0]*1e-3,10.4,color=colors[k],lw=2)
        # Check density
        # ax1.plot(profiles.r*1e-3,premdensity(profiles.r),'k--')
        # ax1.plot(profiles.r*1e-3,profiles.density)
        k+=1

    # Look up AK135
    radius_ak135 = ak135radius()
    vp_ak135 = ak135vp(radius_ak135)
    
    # Plot FVW, PREM and AK135
    ax.plot(profiles.r*1e-3,vp_fvw,color='blue',lw=2,ls=':',label='Ohtaki et al. (2015)')
    ax.vlines(profiles.r[0]*1e-3,vp_fvw[0],10.4,color='blue',lw=2,ls=':')
    ax.plot(profiles.r*1e-3,premvp(profiles.r)*1e-3,color='k',ls='--',label='PREM')
    ax.vlines(profiles.r[0]*1e-3,premvp(profiles.r[0])*1e-3,10.4, 'k', linestyle='--')
    ax.plot(radius_ak135*1e-3,vp_ak135*1e-3,'k',label='ak135')
    ax.vlines(radius_ak135[0]*1e-3,vp_ak135[0]*1e-3,10.4, 'k')

    ax.legend(loc = 0, fontsize=11.5)
    ax.set(xlabel="Radius (km)")
    ax.set(ylabel="P wave speed (km/s)")
    ax.set_xlim([1200,profiles.r.iloc[-1]*1e-3])
    ax.set_ylim([10.25,10.4])
    major_xticks = np.arange(1200,1370,20)
    minor_xticks = np.arange(1200,1370,5)    
    major_yticks = np.arange(10.25,10.4,0.05)
    minor_yticks = np.arange(10.25,10.4,0.01)
    ax.set_xticks(major_xticks)
    ax.set_xticks(minor_xticks, minor=True)
    ax.set_yticks(major_yticks)
    ax.set_yticks(minor_yticks, minor=True)
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)
    ax.tick_params(which='major',length = 7)
    ax.tick_params(which='minor',length = 3.5)

    if saveOn==1:
        saveDir='figures/seismic/'
        if not os.path.exists(saveDir+foldername):
            os.makedirs(saveDir+foldername)
        fig.savefig(saveDir+foldername+"/"+filename+".pdf",format='pdf', dpi=200, bbox_inches='tight')
        fig.savefig(saveDir+foldername+"/"+filename+".png",format='png', dpi=200, bbox_inches='tight')
        print('Figure saved as {}'.format(saveDir+foldername+"/"+filename+".pdf"))
    plt.show()

    return profiles.r, vp_slurry, profiles.density, vp_fvw

# %%
def plot_seismic_dark(layer_thickness, thermal_conductivity,
                             icb_heatflux, csb_heatflux,saveOn,figAspect=0.75):
    w, h = plt.figaspect(figAspect)
    fig, ax = plt.subplots(1,1,figsize=(1.25*w,1.25*h))

    n = layer_thickness.size*thermal_conductivity.size \
        *icb_heatflux.size*csb_heatflux.size

    if n!=1:
        fig_label = [r'high $Le$', r'low $Le$']
        start = 0.2
        stop = 0.5
        cm_subsection = np.linspace(start, stop, n)
        colors = [ cm.gray_r(x) for x in cm_subsection ]

    # Load data
    k=0
    for w,x,y,z in [(w,x,y,z) for w in layer_thickness for x in icb_heatflux for y in csb_heatflux for z in thermal_conductivity]:
        foldername, filename = get_outputDir(w,x,y,z)
        inputDir = "results/{}/{}/".format(foldername,filename)

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

        max_diff = np.max((vp_fvw-vp_slurry*1e-3)/vp_fvw*100)
        print('Maximum difference with Ohtaki et al. (2015) is {:.2f}%'.format(max_diff))
        max_diff = np.max((premvp(profiles.r)-vp_slurry)/premvp(profiles.r)*100)
        print('Maximum difference with PREM is {:.2f}%'.format(max_diff))
        print('Density on slurry side of ICB is {:.2f}'.format(profiles.density[0]))
        density_jump = profiles.density[0] - profiles.density.iloc[-1]
        print('Density jump is {:.2f}'.format(density_jump))
        rho_bod = density_solidFe - profiles.density[0]
        print('Delta rho bod is {:.2f}'.format(rho_bod))
        rho_mod = rho_bod + density_jump
        print('Delta rho mod is {:.2f}'.format(rho_mod))

        # Plot P wave speed
        if n==1:
            ax.plot(profiles.r*1e-3,vp_slurry*1e-3,color='darkgrey',lw=2,label='slurry') #(km/s)
            ax.vlines(profiles.r[0]*1e-3,vp_slurry[0]*1e-3,10.4,color='darkgrey',lw=2)
        else:
            ax.plot(profiles.r*1e-3,vp_slurry*1e-3,color=colors[k],lw=2,label=fig_label[k]) #(km/s)
            ax.vlines(profiles.r[0]*1e-3,vp_slurry[0]*1e-3,10.4,color=colors[k],lw=2)
        # Check density
        # ax1.plot(profiles.r*1e-3,premdensity(profiles.r),'k--')
        # ax1.plot(profiles.r*1e-3,profiles.density)
        k+=1

    # Look up AK135
    radius_ak135 = ak135radius()
    vp_ak135 = ak135vp(radius_ak135)

    ax.plot(profiles.r*1e-3,vp_fvw,color='blue',lw=2,ls=':',label='Ohtaki et al. (2015)')
    ax.vlines(profiles.r[0]*1e-3,vp_fvw[0],10.4,color='blue',lw=2,ls=':')
    ax.plot(profiles.r*1e-3,premvp(profiles.r)*1e-3,color='white',ls='--',label='PREM')
    ax.vlines(profiles.r[0]*1e-3,premvp(profiles.r[0])*1e-3,10.4, 'white', linestyle='--')
    ax.plot(radius_ak135*1e-3,vp_ak135*1e-3,'white',label='ak135')
    ax.vlines(radius_ak135[0]*1e-3,vp_ak135[0]*1e-3,10.4, 'white')

    ax.legend(fontsize=11.5)
    ax.set(xlabel="Radius (km)")
    ax.set(ylabel="P wave speed (km/s)")
    ax.set_xlim([1200,profiles.r.iloc[-1]*1e-3])
    ax.set_ylim([10.1,10.4])
    plt.yticks(np.arange(10.1,10.41,0.1))

    if saveOn==1:
        saveDir='figures/seismic/'
        if not os.path.exists(saveDir+foldername):
            os.makedirs(saveDir+foldername)
        fig.savefig(saveDir+foldername+"/"+filename+"_dark.pdf",format='pdf', dpi=200, bbox_inches='tight')
        fig.savefig(saveDir+foldername+"/"+filename+"_dark.png",format='png', dpi=200, bbox_inches='tight')
        print('Figure saved as {}'.format(saveDir+foldername+"/"+filename+".pdf"))
    plt.show()

    return profiles.r, vp_slurry, profiles.density

# %%
def plot_compare(layer_thickness,csb_heatflux,icb_heatflux,thermal_conductivity,
                 saveOn,saveTag,mol_conc_oxygen_bulk=8.,mol_conc_SSi=8.,
                 self_diffusion=0.98e-8,aspectRatio=0.75,tempMin=5400,
                 tempMax = 5800,xiMin = 6,xiMax = 8,jMin = -3.5e-7,jMax = 0,
                 denMin=11900,denMax = 12250):
    
    w,h= plt.figaspect(aspectRatio)*2
    fig=plt.figure()
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,figsize=(w,h))
    colors = plt.cm.tab10(np.linspace(0,1,layer_thickness.size)) 
    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab10.colors)
    
    for i in (range(layer_thickness.size)):
        foldername,filename = get_outputDir(layer_thickness[i],icb_heatflux,
                                            csb_heatflux,thermal_conductivity)
        inputDir=foldername+"/"+filename
        data_in,data_out,data_profiles=readdata(inputDir)
   
        radius=(data_profiles['r'])*1e-3
        oxygen=data_profiles['oxygen']
        (mass_conc_O,acore) =getcsbmassoxygen(data_in.oxygen_bulk)
        acore=float(acore)
        ax1.plot(radius,data_profiles['temp'],label='_nolegend_')#, color=colors[i])
        ax2.plot(radius,oxygen*acore/aO*100)#, color=colors[i])
        ax3.plot(radius,data_profiles['solidflux'])#,color=colors[i])
        ax4.plot(radius,data_profiles['density'],
                 label='{:.0f} km'.format(layer_thickness[i]*1e-3))#,color=colors[i])
    
    # Liquidus        
    radius_liquidus=np.linspace(icb_radius,icb_radius+400e3)
    temp_liquidus=liquidus(radius_liquidus)
    ax1.plot(radius_liquidus*1e-3,temp_liquidus,'k--',label = 'Liquidus (Davies et al. 2015)')
     
    # PREM
    radius_prem=np.linspace(icb_radius,icb_radius+400e3)
    density_prem=premdensity(radius_prem)
    ax4.plot(radius_prem*1e-3,density_prem,'k--',label='PREM')
    ax4.set(xlabel="Radius (km)",ylabel="Density ($\mathrm{kg m^{-3}}$)")
    
    # Axis titles
    ax1.set(ylabel="Temperature (K)")
    ax2.set(ylabel="Oxygen (mol.%)")
    ax3.set(xlabel="Radius (km)",ylabel="Solid flux ($\mathrm{kg m^{-2} s^{-1}}$)")
    
    # Legend
    ax1.legend(fontsize=11.5,loc=1)
    ax4.legend(fontsize=11.5,loc=1)
    
    # Axis limits
    ax1.set_ylim([tempMin,tempMax])
    ax2.set_ylim([xiMin,xiMax])
    ax3.set_ylim([jMin,jMax])
    ax4.set_ylim([denMin,denMax])
    
    ax1.set_xlim([1220,(icb_radius+400e3)*1e-3])
    
    # Subfigure labels
    ax1.text(1225,tempMax-23,'(a)',fontsize=14)
    ax2.text(1225,xiMax - 0.12,'(b)',fontsize=14)
    ax3.text(1225,jMax - .2e-7,'(c)',fontsize=14)
    ax4.text(1225,denMax - 20,'(d)',fontsize=14)
    
    plt.tight_layout()
    
    if saveOn==1:
        if not os.path.exists('figures/profiles'):
            os.makedirs('figures/profiles')
        saveName=foldername+"_"+filename+saveTag
        plt.savefig('figures/profiles/'+saveName+'.pdf',format='pdf',dpi=200)
        
    plt.show()

# %%
def plot_CD(layer_thickness, thermal_conductivity,
                             icb_heatflux, csb_heatflux,saveOn,figAspect=0.75):
    w, h = plt.figaspect(figAspect)
    # fig, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(1.25*w,1.25*h),sharex=True)
    fig = plt.figure(constrained_layout=True,figsize=(1.25*w,2*h))
    spec = gridspec.GridSpec(ncols=1, nrows=3, figure=fig)
    ax1 = fig.add_subplot(spec[0, 0])
    ax2 = ax1.twinx()
    ax3 = fig.add_subplot(spec[1, 0],sharex=ax1)

    n = layer_thickness.size*thermal_conductivity.size \
        *icb_heatflux.size*csb_heatflux.size

    if n!=1:
        fig_label = [r'high $Le$', r'low $Le$']
        start = 0.2
        stop = 0.5
        cm_subsection = np.linspace(start, stop, n)
        colors = [ cm.gray_r(x) for x in cm_subsection ]

    # Load data
    k=0
    for w,x,y,z in [(w,x,y,z) for w in layer_thickness for x in icb_heatflux for y in csb_heatflux for z in thermal_conductivity]:
        foldername, filename = get_outputDir(w,x,y,z)
        inputDir = "results/{}/{}/".format(foldername,filename)

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

        max_diff = np.max((vp_fvw-vp_slurry*1e-3)/vp_fvw*100)
        print('Maximum difference with Ohtaki et al. (2015) is {:.2f}%'.format(max_diff))
        max_diff = np.max((premvp(profiles.r)-vp_slurry)/premvp(profiles.r)*100)
        print('Maximum difference with PREM is {:.2f}%'.format(max_diff))
        print('Density on slurry side of ICB is {:.2f}'.format(profiles.density[0]))
        density_jump = profiles.density[0] - profiles.density.iloc[-1]
        print('Density jump is {:.2f}'.format(density_jump))
        rho_bod = density_solidFe - profiles.density[0]
        print('Delta rho bod is {:.2f}'.format(rho_bod))
        rho_mod = rho_bod + density_jump
        print('Delta rho mod is {:.2f}'.format(rho_mod))

        # Plot P wave speed
        if n==1:
            ax3.plot(profiles.r*1e-3,vp_slurry*1e-3,color='darkgrey',lw=2,label='slurry density') #(km/s)
            ax3.vlines(profiles.r[0]*1e-3,vp_slurry[0]*1e-3,10.4,color='darkgrey',lw=2)
        else:
            ax3.plot(profiles.r*1e-3,vp_slurry*1e-3,color=colors[k],lw=2,label=fig_label[k]) #(km/s)
            ax3.vlines(profiles.r[0]*1e-3,vp_slurry[0]*1e-3,10.4,color=colors[k],lw=2)
        # Check density
        # ax3.plot(profiles.r*1e-3,premdensity(profiles.r),'k--')
        # ax3.plot(profiles.r*1e-3,profiles.density)
        k+=1
    
    # Plot temperature and composition
    ax1.plot(profiles.r*1e-3,profiles.temp,lw=2,color='red',label='temperature')
    (mass_conc_O,acore) = getcsbmassoxygen(inputs.oxygen_bulk)
    acore=float(acore)
    ax2.plot(profiles.r*1e-3,profiles.oxygen*acore/aO*100,lw=2,color='blue',label='oxygen')    
    
    # Plot FVW, PREM and AK135
    ax3.plot(profiles.r*1e-3,vp_fvw,color='black',lw=2,ls=':',label='Ohtaki et al. (2015)')
    ax3.vlines(profiles.r[0]*1e-3,vp_fvw[0],10.4,color='black',ls=':',lw=2)
    ax3.plot(profiles.r*1e-3,premvp(profiles.r)*1e-3,color='k',ls='--',label='PREM')
    ax3.vlines(profiles.r[0]*1e-3,premvp(profiles.r[0])*1e-3,10.4, 'k', linestyle='--')
    # ax3.plot(radius_ak135*1e-3,vp_ak135*1e-3,'k',label='ak135')
    # ax3.vlines(radius_ak135[0]*1e-3,vp_ak135[0]*1e-3,10.4, 'k')

    ax3.legend(fontsize=11.5, loc='upper right')
    ax1.set(ylabel="Temperature (K)")
    ax1.spines["right"].set_edgecolor('red')
    ax1.spines["right"].set_edgecolor('blue')
    ax1.yaxis.label.set_color('red')
    ax2.set(ylabel="Oxygen (mol.%)")
    ax2.spines["left"].set_edgecolor('red')
    ax2.spines["right"].set_edgecolor('blue')
    ax2.yaxis.label.set_color('blue')
    ax2.tick_params(axis='y', colors='blue')
    ax1.tick_params(axis='y', which='both', colors='red')
    ax3.set(xlabel="Radius (km)")
    ax3.set(ylabel="P wave speed (km/s)")
    ax3.set_xlim([1200,profiles.r.iloc[-1]*1e-3])
    ax3.set_ylim([10.25,10.4])
    major_xticks = np.arange(1200,1370,20)
    minor_xticks = np.arange(1200,1370,5)    
    ax1.set_xticks(major_xticks)
    ax1.set_xticks(minor_xticks, minor=True)
    ax2.set_xticks(major_xticks)
    ax2.set_xticks(minor_xticks, minor=True)    
    ax3.set_xticks(major_xticks)
    ax3.set_xticks(minor_xticks, minor=True)
    major_yticks = np.arange(5500,5751,100)
    minor_yticks = np.arange(5500,5751,20)
    ax1.set_yticks(major_yticks)
    ax1.set_yticks(minor_yticks, minor=True)
    ax1.grid(which='minor', alpha=0.2)
    ax1.grid(which='major', alpha=0.5)
    ax1.tick_params(which='major',length = 7)
    ax1.tick_params(which='minor',length = 3.5)    
    # major_yticks = np.arange(6.5,8.1,0.5)
    # minor_yticks = np.arange(6.5,8.1,0.1)
    # ax2.set_yticks(major_yticks)
    # ax2.set_yticks(minor_yticks, minor=True)
    # ax2.grid(which='minor', alpha=0.2)
    # ax2.grid(which='major', alpha=0.5)
    # ax2.tick_params(which='major',length = 7)
    # ax2.tick_params(which='minor',length = 3.5)
    major_yticks = np.arange(10.25,10.4,0.05)
    minor_yticks = np.arange(10.25,10.4,0.01)
    ax3.set_yticks(major_yticks)
    ax3.set_yticks(minor_yticks, minor=True)
    ax3.grid(which='minor', alpha=0.2)
    ax3.grid(which='major', alpha=0.5)
    ax3.tick_params(which='major',length = 7)
    ax3.tick_params(which='minor',length = 3.5)

    if saveOn==1:
        saveDir='figures/seismic/'
        if not os.path.exists(saveDir+foldername):
            os.makedirs(saveDir+foldername)
        fig.savefig(saveDir+foldername+"/"+filename+"_CD.pdf",format='pdf', dpi=200, bbox_inches='tight')
        fig.savefig(saveDir+foldername+"/"+filename+"_CD.png",format='png', dpi=200, bbox_inches='tight')
        print('Figure saved as {}'.format(saveDir+foldername+"/"+filename+"_CD.pdf"))
    plt.show()

    return profiles.r, vp_slurry, profiles.density, vp_fvw