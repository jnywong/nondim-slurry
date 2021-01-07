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

#def getparameters(layer_thickness,csb_heatflux,thermal_conductivity):
#                  icb_heatflux,ic_age):        
import numpy as np
from scipy.constants import gas_constant
from scipy import integrate, interpolate

from slurpy.coreproperties import icb_radius,density_liquidO, \
density_solidFe,heat_capacity,latent_heat,aO,aFe,aSSi,deltaV_solidFe_liquidFe, \
gigayear, alphaT, alphaXi, bulk_modulus, cmb_radius, gruneisen
from slurpy.lookup import liquidus, premgravity, premdensity, premvp, ohtaki

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
def getcsbheatflux(St,freezing_speed,density0,csb_radius):
    csb_heatflux = St*density0*freezing_speed*latent_heat*4*np.pi*csb_radius**2*1e-12
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

# IC age given ICB speed
def geticage(icb_speed):
    ic_age=icb_radius/(3*icb_speed*gigayear)
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
def getLip(csb_radius,model='prem'):
    gravity=premgravity(csb_radius)
    if model == 'prem':
        density=premdensity(csb_radius)
    elif model == 'ohtaki':
        _,density = ohtaki(csb_radius)
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
def getLewis(thermal_conductivity,self_diffusion,density0):
    Le=thermal_conductivity/(density0*
                                heat_capacity*self_diffusion)
    return Le

# Peclet number
def getPeclet(freezing_speed,csb_radius,self_diffusion):
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
    csb_heatflux=getcsbheatflux(St,freezing_speed,density0,csb_radius)
    thermal_conductivity=getthermalconductivity(Le, density0)
    return (icb_heatflux,csb_heatflux,thermal_conductivity)

# %% Postprocessing
def slurrydensity(radius,temp,xi,solidflux,mol_conc_oxygen_bulk, \
                  sedimentation_constant,mol_conc_SSi=8,model='prem'):
    # Density profile across slurry layer    
    csb_radius=radius[-1]
    if model=='prem':
        csb_density=premdensity(csb_radius)
    elif model =='ohtaki':
        _,csb_density=ohtaki(csb_radius)
    # Solid fraction
    Kphi=getKphi(sedimentation_constant,radius,mol_conc_oxygen_bulk)
    phi=(-solidflux/Kphi)**(3/5)
    # Density fluctuations
    deltaV_solid_liquid=getchangevolmeltingFe(csb_density)
    temp_denFluc=-csb_density*alphaT*(temp-temp[-1])
    xi_denFluc=-csb_density*alphaXi*(xi-xi[-1])
    phi_denFluc=csb_density*(csb_density*deltaV_solid_liquid+alphaXi*xi)*(phi-phi[-1])
    density_fluc= temp_denFluc+xi_denFluc+phi_denFluc
    # Hydrostatic density
    slurry_gravity=premgravity(radius)
    density_hydrostatic=1/(slurry_gravity*(radius-csb_radius)/bulk_modulus+ \
                           1/csb_density) 
    # Total density
    density_slurry=density_hydrostatic+density_fluc
    return (density_slurry, phi, temp_denFluc, xi_denFluc, phi_denFluc,density_fluc)

# %%    
def adiabat(oc_radius,csb_temp,n):
    # Construct adiabat across outer core given CSB temperature
    # Gravity
    oc_gravity=premgravity(oc_radius)
    # Seismic parameter
    seismic_parameter=premvp(oc_radius)**2
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
    temp_ad1 = adiabat(radius_oc1,liquidus(csb_radius1),50)
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
    gravity_psi = premgravity(radius_psi)
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
    density_oc = premdensity(radius_oc)
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