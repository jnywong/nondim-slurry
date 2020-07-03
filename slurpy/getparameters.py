#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 11:19:34 2018

@author: jennywong
"""

###############################################################################
# GETPARAMETERS.PY                                                            #  
###############################################################################

# Description: reference parameters to be passed into the BVP for the
#               spherical steady state slurry
    
import numpy as np
from scipy.constants import gas_constant

from slurpy.coreproperties import icb_radius,density_liquidO, \
density_solidFe,heat_capacity,latent_heat,aO,aFe,aSSi,deltaV_solidFe_liquidFe, \
gigayear, alphaT, alphaXi, bulk_modulus, cmb_radius, gruneisen
from slurpy.lookup import liquidus, premgravity, premdensity, premvp

# CSB radius
def getcsbradius(layer_thickness):
    csb_radius=icb_radius+layer_thickness
    return csb_radius

# CSB temperature
def getcsbtemp(layer_thickness):
    csb_radius= getcsbradius(layer_thickness)
    csb_temp=liquidus(csb_radius)
    return csb_temp

# CSB temperature gradient
def getcsbtempgrad(layer_thickness,csb_heatflux,thermal_conductivity):
    csb_radius=getcsbradius(layer_thickness)
    csb_temp_grad=csb_heatflux*1e12/(-4*thermal_conductivity*np.pi*csb_radius**2)
    return csb_temp_grad

def getcsbmassoxygen(mol_conc_oxygen_bulk, mol_conc_SSi=8):   
    mol_conc_O=mol_conc_oxygen_bulk*1e-2
    mol_conc_SSi = mol_conc_SSi*1e-2
    mol_conc_Fe = 1-mol_conc_O-mol_conc_SSi
    acore = mol_conc_Fe*aFe+mol_conc_SSi*aSSi+mol_conc_O*aO # atomic weight of core material
    # CSB oxygen concentration
    mass_conc_O = mol_conc_O*aO/acore
    return mass_conc_O,acore
         
# Change in specific volume upon melting
def getchangevolmelting(mol_conc_oxygen_bulk,density0,mol_conc_SSi=8):
    mass_conc_O = getcsbmassoxygen(mol_conc_oxygen_bulk, mol_conc_SSi)[0]
    deltaV_liquidFeO_solidFe = (1-mass_conc_O)/density0 + \
    mass_conc_O/density_liquidO - 1/density_solidFe
    return deltaV_liquidFeO_solidFe

# Change in specific volume of liquid and solid iron
def getchangevolmeltingFe(density0):
    # Change in specific volume between solid and liquid
    deltaV_solid_liquid = 1/density0 - 1/density_solidFe
    return deltaV_solid_liquid

# Freezing speed
def getfreezingspeed(icb_heatflux):
    freezing_speed=icb_heatflux*1e12/(4*np.pi*icb_radius**2*density_solidFe*latent_heat)
    return freezing_speed

# ICB heat flux
def geticbheatflux(freezing_speed):
    icb_heatflux = freezing_speed*(4*np.pi*icb_radius**2*density_solidFe*latent_heat)*1e-12
    return icb_heatflux

# CSB heat flux
def getcsbheatflux(St,freezing_speed,density0,csb_radius,icb_heatflux):
    # csb_heatflux = St*density0*freezing_speed*latent_heat*4*np.pi*csb_radius**2*1e-12
    csb_heatflux = St*icb_heatflux*csb_radius**2/icb_radius**2
    return csb_heatflux

# Thermal conductivity
def getthermalconductivity(Le,density0,self_diffusion=0.98e-8):
    thermal_conductivity=Le*density0*heat_capacity*self_diffusion
    return thermal_conductivity

# Initial guess for snow speed
def getinitialsnowspeed(ic_age):
    initial_snow_speed = icb_radius/(3*ic_age*gigayear) # assume volumetric growth rather than linear growth
    return initial_snow_speed

# Core cooling rate
def getcorecoolingrate(core_cooling_rate):
    core_cooling_rate_out=core_cooling_rate/gigayear
    return core_cooling_rate_out

# Snow speed given IC age
def getsnowspeed(ic_age):
    snow_speed=icb_radius/(3*ic_age*gigayear)
    return snow_speed

# IC age given snow speed
def geticage(snow_speed):
    ic_age=icb_radius/(3*snow_speed*gigayear)
    return ic_age

# Freezing speed given Pe
def getfreezingspeedfromPe(Pe,csb_radius,self_diffusion):
    freezing_speed = Pe*self_diffusion/csb_radius
    return freezing_speed

#%%
# Kphi - prefactor in solid flux eq (2.9.1g)
def getKphi(sedimentation_constant, radius, mol_conc_oxygen_bulk, mol_conc_SSi=8):
    gravity=premgravity(radius)
    density=premdensity(radius)
    deltaV_solid_liquid = getchangevolmeltingFe(density[-1])
    Kphi=sedimentation_constant*gravity*density*deltaV_solid_liquid
    return Kphi

def getphi(Kphi,solidflux):
    phi = (-solidflux/Kphi)**(3/5)
    return phi

#%% Scalings
# Temp
def get_tempScale(csb_heatflux,csb_radius,density0,freezing_speed):
    scale_temp = csb_heatflux*1e12/(4*np.pi*csb_radius**2*density0* \
                                    heat_capacity*freezing_speed)
    return scale_temp

# Solid flux
def get_jScale(freezing_speed):
    scale_j = freezing_speed*density_solidFe
    return scale_j

#%% Dimensionless parameters
# Lip
def getLip(csb_radius):
    gravity=premgravity(csb_radius)
    density=premdensity(csb_radius)
    Lip=deltaV_solidFe_liquidFe*gravity*density*csb_radius/latent_heat
    return Lip,gravity,density

# Lix
def getLix(mass_conc_O):
    L2=1000*gas_constant*mass_conc_O/(aO*heat_capacity)
    return L2

# Stefan number
def getStefan(icb_heatflux,csb_heatflux,csb_radius):
    St=csb_heatflux*icb_radius**2/(icb_heatflux*csb_radius**2)
    return St

# Lewis number
def getLewis(thermal_conductivity,density0,self_diffusion=0.98e-8):
    Le=thermal_conductivity/(density0*
                                heat_capacity*self_diffusion)
    return Le

# Peclet number
def getPeclet(freezing_speed,csb_radius,self_diffusion=0.98e-8):
    Pe=freezing_speed*csb_radius/(self_diffusion)
    return Pe

#%% Dimensional parameters
def getdimensional(layer_thickness,Pe,St,Le, \
                   mol_conc_oxygen_bulk=8., mol_conc_SSi=8., \
                       self_diffusion=0.98e-8):
    csb_radius = getcsbradius(layer_thickness)
    density0=premdensity(csb_radius)
    mass_conc_O,acore=getcsbmassoxygen(mol_conc_oxygen_bulk, mol_conc_SSi)
    freezing_speed = getfreezingspeedfromPe(Pe,csb_radius,self_diffusion)
    icb_heatflux=geticbheatflux(freezing_speed)
    csb_heatflux=getcsbheatflux(St,freezing_speed,density0,csb_radius,icb_heatflux)
    thermal_conductivity=getthermalconductivity(Le, density0)
    return (icb_heatflux,csb_heatflux,thermal_conductivity)

#%% Directory name
def getdirectory(layer_thickness,icb_heatflux,csb_heatflux,thermal_conductivity, \
                 mol_conc_oxygen_bulk=8.,self_diffusion=0.98e-8):
    csb_radius = getcsbradius(layer_thickness)
    mass_conc_O,acore = getcsbmassoxygen(mol_conc_oxygen_bulk)
    freezing_speed = getfreezingspeed(icb_heatflux)
    Lip,_,density0 = getLip(csb_radius)
    Lix = getLix(mass_conc_O)
    Le = getLewis(thermal_conductivity,self_diffusion,density0)
    Pe = getPeclet(freezing_speed,csb_radius,self_diffusion)
    St = getStefan(icb_heatflux,csb_heatflux,csb_radius)
    str1=str(np.round(Le,2)).replace('.','_')
    str2=str(np.round(Lip,2)).replace('.','_')
    str3=str(np.round(Lix,3)).replace('.','_')
    str4=str(np.round(Pe,2)).replace('.','_')
    str5=str(np.round(St,2)).replace('.','_')
    directory="results/Le_{}/Lip_{}_Lix_{}_Pe_{}_St_{}/".format(str1,str2,str3,str4,str5)     
    
    return (directory,Lip,Lix,Le,Pe,St)

# Critical St 
def getcriticalSt(layer_thickness):
    St_crit = icb_radius**2/(icb_radius+layer_thickness)**2
    
    return St_crit
