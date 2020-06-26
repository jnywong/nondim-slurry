#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 11:08:06 2020

@author: jennywong
"""

###############################################################################
# POSTPROCESS.PY                                                              #  
###############################################################################

# Description: routines to postprocess results from the spherical steady 
#               state slurry

import numpy as np
import pandas as pd
from scipy import integrate, interpolate

import slurpy.getparameters as gp
import slurpy.lookup as lp
from slurpy.coreproperties import icb_radius, density_solidFe,heat_capacity, \
    latent_heat, gigayear, alphaT, alphaXi, bulk_modulus, cmb_radius, gruneisen

# %% Postprocessing
def slurrydensity(radius,temp,xi,solidflux,mol_conc_oxygen_bulk, \
                  sedimentation_constant,mol_conc_SSi=8):
    # Density profile across slurry layer    
    csb_radius=radius[-1]
    csb_density=lp.premdensity(csb_radius)
    # Solid fraction
    Kphi=gp.getKphi(sedimentation_constant,radius,mol_conc_oxygen_bulk)
    phi=(-solidflux/Kphi)**(3/5)
    # Density fluctuations
    deltaV_solid_liquid=gp.getchangevolmeltingFe(csb_density)
    temp_denFluc=-csb_density*alphaT*(temp-temp[-1])
    xi_denFluc=-csb_density*alphaXi*(xi-xi[-1])
    phi_denFluc=csb_density*(csb_density*deltaV_solid_liquid+alphaXi*xi)*(phi-phi[-1])
    density_fluc= temp_denFluc+xi_denFluc+phi_denFluc
    # Hydrostatic density
    slurry_gravity=lp.premgravity(radius)
    density_hydrostatic=1/(slurry_gravity*(radius-csb_radius)/bulk_modulus+ \
                           1/csb_density) 
    # Total density
    density_slurry=density_hydrostatic+density_fluc
    return (density_slurry, phi, temp_denFluc, xi_denFluc, phi_denFluc,density_fluc)

# %%    
def adiabat(oc_radius,csb_temp,n):
    # Construct adiabat across outer core given CSB temperature
    # Gravity
    oc_gravity=lp.premgravity(oc_radius)
    # Seismic parameter
    seismic_parameter=lp.premvp(oc_radius)**2
    # Pre-allocate 
    temp_adiabat=np.zeros(n)
    # Integrand for the adiabat
    integrand=oc_gravity*gruneisen/seismic_parameter
    # Boundary condition
    temp_adiabat[0] = csb_temp
    # Integrate
    for i in range(1,n):
        temp_adiabat[i]=csb_temp*np.exp(-integrate.simps(integrand[0:i+1], \
                    oc_radius[0:i+1]))  
    return temp_adiabat

#%% Calculate difference in CMB temp after 1Ga by constructing adiabats
# anchored at present day and after 1Ga
def get_cooling(icb_speed,csb_temp, radius_oc, csb_radius, cmb_radius):
    delta_t = gigayear
    
    temp_ad0 = adiabat(radius_oc,csb_temp,50)
    cmb_temp0 = temp_ad0[-1]
    print('CMB temp present day is {:.0f}K'.format(cmb_temp0))
    
    csb_radius1 = csb_radius+icb_speed*delta_t
    radius_oc1 = np.linspace(csb_radius1,cmb_radius)
    temp_ad1 = adiabat(radius_oc1,lp.liquidus(csb_radius1),50)
    cmb_temp1 = temp_ad1[-1]
    print('CMB temp after 1Ga is {:.0f}K'.format(cmb_temp1))
    cooling_rate = (cmb_temp1 - cmb_temp0)/delta_t
    cooling_rate_Ga = cooling_rate*gigayear
    print('Cooling rate is {:.2f}K/Ga'.format(cooling_rate_Ga))
    return  cooling_rate, cmb_temp0, temp_ad0

#%% Heat flux across core-mantle boundary
def heatflux(radius,temp,xi,solidflux,phi,temp_grad,xi_grad,density_slurry, \
             icb_speed,icb_heatflux,layer_thickness,thermal_conductivity, \
             csb_heatflux,n):

    csb_radius = radius[-1]      

    # GRAVITATIONAL POWER
    # Gravitational potential
    radius_psi = np.linspace(icb_radius,cmb_radius)
    gravity_psi = lp.premgravity(radius_psi)
    psi = np.zeros(radius_psi.size)
    for i in range(radius_psi.size-1):
        psi[i] = -integrate.simps(gravity_psi[i:],radius_psi[i:])  

    # Slurry
    f = interpolate.interp1d(radius_psi,psi)
    psi_slurry = f(radius)    
    oxygen_mass_slurry = 4*np.pi*density_solidFe*layer_thickness**2*icb_speed*xi[-1]
    mass_slurry = integrate.simps(density_slurry*4*np.pi*radius**2,radius)
    Qg_slurry = -alphaXi*oxygen_mass_slurry/mass_slurry* \
        integrate.simps(density_slurry*psi_slurry*4*np.pi*radius**2,radius)
#    print('Change in oxygen mass (slurry) is {}'.format(oxygen_mass_slurry))
#    print('Mass of slurry is {}'.format(mass_slurry))
    print('Qg slurry is {:.2f}TW'.format(Qg_slurry*1e-12))

    # Outer core
    radius_oc = np.linspace(csb_radius,cmb_radius)
    psi_oc = f(radius_oc)
    surf = -psi[0]*alphaXi*density_slurry[-1]*xi[-1]*icb_speed*4*np.pi*csb_radius**2

    # Mass of OC
    density_oc = lp.premdensity(radius_oc)
    mass_oc = integrate.simps(density_oc*4*np.pi*radius_oc**2,radius_oc)
    oxygen_mass_oc = - density_slurry[-1]*icb_speed*4*np.pi*csb_radius**2*xi[-1]/mass_oc
    bulk = integrate.simps(alphaXi*psi_oc*oxygen_mass_oc*4*np.pi*radius_oc**2,radius_oc)
#    print('Change in oxygen mass (outer core) is {}'.format(oxygen_mass_oc))
#    print('Mass of outer core is {}'.format(mass_oc))
    print('Qg surface term in outer core is {:.2f}TW'.format(surf*1e-12))
    print('Qg bulk term in outer core is {:.5f}TW'.format(bulk*1e-12))    
    
    # Total gravitational energy
    Qg_oc = surf+bulk
    Qg = Qg_slurry + Qg_oc
    print('Total Qg is {:.2f}TW'.format(Qg*1e-12))    
    
    # LATENT HEAT
    Ql=4*np.pi*icb_radius**2*density_solidFe*icb_speed*latent_heat
    print('Total Ql is {:.2f}TW'.format(Ql*1e-12))  
    
    # SECULAR COOLING    
    # Cooling rate
    csb_temp = temp[-1]                    
    cooling_rate,cmb_temp,temp_ad = get_cooling(icb_speed,csb_temp,radius_oc,
                                                    csb_radius, cmb_radius)
    
    # Outer core
    Qs_oc = - heat_capacity*cooling_rate/cmb_temp \
            *integrate.simps(density_oc*temp_ad*4*np.pi*radius_oc**2,radius_oc)
    print('Qs in outer core is {:.2f}TW'.format(Qs_oc*1e-12))
    
    # Slurry
    Qs_slurry = csb_heatflux - Qg_slurry - Ql
    print('Qs in slurry {:.2f}TW'.format(Qs_slurry*1e-12))
    Qs = Qs_slurry+Qs_oc
    print('Total Qs is {:.2f}TW'.format(Qs*1e-12))
    
    # CMB heat flux
    Q_cmb = csb_heatflux+Qs_oc+Qg_oc
    print('Qcmb is is {:.2f}TW'.format((Q_cmb)*1e-12))
    
    return (Q_cmb, Qs, Qs_slurry, Qs_oc, Ql, Qg, Qg_oc, Qg_slurry, cooling_rate,
            cmb_temp, temp_ad)

#%% Stable layer
def stablelayer(density_fluc,radius):
    density_fluc_grad = np.gradient(density_fluc,radius)
    idx = np.argwhere(density_fluc_grad>0)
    if idx.size!=0:
        # Unstable or partiall stable?
        unstable_range = radius[idx[-1][0]] - radius[idx[0][0]]
        if unstable_range > 100e3:
            stable = 0 # unstable
        else:
            stable = 2 # partially stable
    else:
        stable = 1 # stable
        
    return stable

#%% Aggregate data for regime diagram
def regime(layer_thickness,thermal_conductivity,csb_heatfluxes,
                icb_heatfluxes,mol_conc_oxygen_bulk=8.,self_diffusion=0.98e-8):
    # Directory
    nPe = icb_heatfluxes.size
    nSt = csb_heatfluxes.size
    Pe = np.zeros((nPe, nSt))
    St = np.zeros((nPe, nSt))
    density_jump = np.zeros((nPe, nSt))
    Q_cmb = np.zeros((nPe, nSt))
    stable = np.zeros((nPe, nSt))
    St_crit = gp.getcriticalSt(layer_thickness)
    for i in range(nPe):
        for j in range(nSt):
            (inputDir,_,_,_,outPe,outSt) = gp.getdirectory(layer_thickness,icb_heatfluxes[i], 
                                       csb_heatfluxes[j], thermal_conductivity)
            try:
                data_input=pd.read_csv(inputDir+"inputs.csv",index_col=False)           
            except FileNotFoundError:
                Pe[i,j] = outPe
                St[i,j] = outSt
                density_jump[i,j] = np.nan
                Q_cmb[i,j] = np.nan
                stable[i,j] = np.nan
                if outSt>St_crit:
                    print('No solution: {}'.format(inputDir))
                continue
            
            data_output=pd.read_csv(inputDir+"outputs.csv",index_col=False)
            Pe[i,j] = np.float(data_input['Pe'])
            St[i,j] = np.float(data_input['St'])
            density_jump[i,j] = np.float(data_output['density_jump'])
            Q_cmb[i,j] = np.float(data_output['Q_cmb'])
            stable[i,j] = np.float(data_output['stable'])     
            
    return (Pe, St, density_jump, Q_cmb, stable)            
    