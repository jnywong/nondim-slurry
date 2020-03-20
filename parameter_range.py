#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 10:57:34 2019

@author: wong
"""
import numpy as np

from refparameters import getL1, getL2, getStefan, getLewis, getPeclet, \
    getcsbradius, getfreezingspeed, getcsbmassoxygen, getchangevolmelting, \
    geticbheatflux
from coreproperties import deltaV_solidFe_liquidFe, \
    latent_heat
from lookup import premdensity

def getdimensionless(layer_thickness,icb_heatflux,csb_heatflux,
                     thermal_conductivity,mol_conc_oxygen_bulk =8.,
                     mol_conc_SSi = 8.,self_diffusion= 0.98e-8):
    csb_radius=getcsbradius(layer_thickness)
    mass_conc_O,acore=getcsbmassoxygen(mol_conc_oxygen_bulk,mol_conc_SSi)
    density0 = premdensity(csb_radius)
    freezing_speed=getfreezingspeed(icb_heatflux)
    
    L1,_,density0=getL1(csb_radius)
    L2=float(getL2(mass_conc_O))
    St=float(getStefan(icb_heatflux,csb_heatflux,csb_radius))
    Le=float(getLewis(thermal_conductivity,self_diffusion,density0))
    Pe=float(getPeclet(freezing_speed,csb_radius,self_diffusion))
    
    L1 = float(L1)
    
    return (L1,L2,Pe,St,Le)

def getdimension(layer_thickness,thermal_conductivity, Pe, St,
                     mol_conc_oxygen_bulk =8.,
                     mol_conc_SSi = 8.,self_diffusion= 0.98e-8):
    
    csb_radius=getcsbradius(layer_thickness)
    density0 = premdensity(csb_radius)
    mass_conc_O,acore=getcsbmassoxygen(mol_conc_oxygen_bulk,mol_conc_SSi)
    deltaV_liquidFeO_solidFe=getchangevolmelting(mol_conc_oxygen_bulk,mol_conc_SSi,density0)    
    # ICB heat flux (TW)
    freezing_speed = Pe*self_diffusion*deltaV_liquidFeO_solidFe \
                                        /(csb_radius*deltaV_solidFe_liquidFe)
    icb_heatflux = geticbheatflux(freezing_speed)
    # CSB heat flux (TW)
    csb_heatflux = St*4*np.pi*csb_radius**2*density0*freezing_speed \
        *latent_heat*1e-12
    
    return (icb_heatflux,csb_heatflux)