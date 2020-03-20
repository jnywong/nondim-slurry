#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 10:22:50 2018

@author: jennywong
"""

# Contains postprocessing functions:
#   SLURRYDENSITY(Z,TEMP,XI,J) - finds the density profile across the slurry 
#                                   layer
#   ADIABAT(CSB_RADIUS,CSB_TEMP) - construct adiabat across OC given CSB
#                                   temperature
#   HEATFLUX(Z,TEMP,XI,J,PHI,TEMP_GRAD,XI_GRAD,DENSITY_SLURRY,ICB_SPEED,
#   SNOW_SPEED) - CMB heat flux

    
import numpy as np

from refparameters import getcsbradius
from lookup import premdensity,premgravity,liquidus
from coreproperties import alphaT, alphaXi, icb_radius, \
cmb_radius, bulk_modulus, latent_heat,density_solidFe,heat_capacity, \
gigayear
from refparameters import getchangevolmeltingFe, getKphi
from scipy import integrate, interpolate


# %%    
def slurrydensity(radius,temp,xi,solidflux,layer_thickness,mol_conc_oxygen_bulk, \
                  mol_conc_SSi,sedimentation_constant):
    # Density profile across slurry layer    
    csb_radius=getcsbradius(layer_thickness)
    csb_density=premdensity(csb_radius)
    
    # Solid fraction
    Kphi=getKphi(sedimentation_constant,radius,mol_conc_oxygen_bulk,mol_conc_SSi)
    phi=(-solidflux/Kphi)**(3/5)
    
    # Density fluctuations
    deltaV_solid_liquid=getchangevolmeltingFe(csb_density)
    temp_fluc=-csb_density*alphaT*(temp-temp[-1])
    xi_fluc=-csb_density*alphaXi*(xi-xi[-1])
    phi_fluc=csb_density*(csb_density*deltaV_solid_liquid+alphaXi*xi)*(phi-phi[-1])
    density_fluc= temp_fluc+xi_fluc+phi_fluc
    
    # Hydrostatic density
    slurry_gravity=premgravity(radius)
    density_hydrostatic=1/(slurry_gravity*(radius-csb_radius)/bulk_modulus+ \
                           1/csb_density) 
    
    # Total density
    density_slurry=density_hydrostatic+density_fluc
    return (density_slurry, phi, temp_fluc, xi_fluc, phi_fluc,density_fluc)

# %%    
def adiabat(oc_radius,csb_temp,n):
    # Construct adiabat across outer core given CSB temperature
    from coreproperties import gruneisen
#    oc_radius=np.linspace(float(csb_radius),cmb_radius,n)
    # Gravity
    from lookup import premgravity
    oc_gravity=premgravity(oc_radius)
    # Seismic parameter
    from lookup import vpspeed
    seismic_parameter=vpspeed(oc_radius)**2
    # Pre-allocate 
    temp_adiabat=np.zeros(n)
    # Integrand for the adiabat
    integrand=oc_gravity*gruneisen/seismic_parameter
    # Boundary condition
    temp_adiabat[0] = csb_temp
    # Integrate
    from scipy import integrate
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
    return  cooling_rate, cmb_temp0, temp_ad0

#%% Heat flux across core-mantle boundary
def heatflux(z,temp,xi,solidflux,phi,temp_grad,xi_grad,density_slurry, \
             snow_speed,freezing_speed,icb_heatflux,layer_thickness,thermal_conductivity, \
             csb_heatflux,n):

    radius=z
    csb_radius = radius[-1]      
    icb_speed = snow_speed+freezing_speed

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
    cooling_rate_Ga = cooling_rate*gigayear
    print('Cooling rate is {:.2f}K/Ga'.format(cooling_rate_Ga))
    
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
   