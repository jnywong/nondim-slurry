#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 10:03:00 2018

@author: jennywong
"""
import matplotlib.pyplot as plt

def plotprofiles(z,temp,xi,j,density,acore,saveOn,outputDir):

    # Radius
    from coreproperties import icb_radius
    radius=(z+icb_radius)*1e-3
    
    # Temperature
    fig1=plt.figure() # figsize=(3.54,2.66)
    ax1 = plt.subplot()
    ax1.plot(radius,temp)
    plt.xlabel("Radius (km)")
    plt.ylabel("Temperature (K)")
#    plt.tight_layout()
    plt.show()
    
    # Oxygen concentration
    from coreproperties import aO
#    from batch import mol_conc_oxygen_bulk,mol_conc_SSi
    fig2=plt.figure()
    ax2 = plt.subplot()
    ax2.plot(radius,xi*acore/aO*100)
#    plt.ylim((round(xi.min()*acore/aO*100,3),mol_conc_oxygen_bulk))
    plt.ylim((round(xi.min()*acore/aO*100,3),xi.max()*acore/aO*100))
    plt.xlabel("Radius (km)")
    plt.ylabel("Oxygen concentration (mol.%)")
    plt.show()
    
    # Solid flux
    fig3=plt.figure()
    ax3 = plt.subplot()
    ax3.plot(radius,j)
    plt.xlabel("Radius (km)")
    plt.ylabel("Solid flux ($\mathrm{kg m^{-2} s^{-1}}$)")
    ax3.yaxis.major.formatter._useMathText = True # change axis label 1e<x> to 10^<x>
    plt.show()
    
    # Density
    fig4=plt.figure()
    ax4 = plt.subplot()
    from lookup import premdensity
    ax4.plot(radius,density)
    # PREM density
    density_prem=premdensity(radius*1e3)
    ax4.plot(radius,density_prem,'k--')
    plt.xlabel("Radius (km)")
    plt.ylabel("Density ($\mathrm{kg m^{-3}}$)")
    ax4.legend(["slurry","PREM"])
    plt.show()
    
    # Save plots
    if saveOn==1:
        fig1.savefig(outputDir+"temp.eps",format='eps', dpi=1000)
        fig2.savefig(outputDir+"xi.eps",format='eps', dpi=1000)
        fig3.savefig(outputDir+"j.eps",format='eps', dpi=1000)
        fig4.savefig(outputDir+"density.eps",format='eps', dpi=1000)
    