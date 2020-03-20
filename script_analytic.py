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
        'size'   : 16}

matplotlib.rc('font', **font)
#------------------------------------------------------------------------------


#outputDir = 'results/highLe/L10.16_L20.02_Le1180.89_St1.59_Pe1959.06/'
#outputDir = 'results/lowLe/L10.16_L20.02_Le354.27_St1.59_Pe1959.06/'

#outputDir = 'results/highLe/L10.22_L20.02_Le1195.97_St1.13_Pe2316.16/' # d=400km
saveName = 'L10.16_L20.02_Le1180.89_St1.59_Pe1959.06' # d=150km
#saveName = 'L10.22_L20.02_Le1195.97_St1.13_Pe2316.16' # d=400km
    
outputDir = 'results/highLe/{}/'.format(saveName)

inputs,outputs,profiles = read_sphdata(outputDir)

layer_thickness = pd.to_numeric(inputs.layerThickness.iloc[0])
icb_heatflux = pd.to_numeric(inputs.icb_heatflux.iloc[0])
csb_heatflux = pd.to_numeric(inputs.csb_heatflux.iloc[0])
thermal_conductivity = pd.to_numeric(inputs.thermalConductivity.iloc[0])

csb_radius = icb_radius + layer_thickness
r_out = (profiles.z)/csb_radius

freezing_speed=getfreezingspeed(icb_heatflux)
F = pd.to_numeric(outputs.F.iloc[0])
snow_speed = pd.to_numeric(outputs.snowSpeed.iloc[0])
icb_speed = snow_speed + freezing_speed

(L1,L2,Pe,St,Le) = getdimensionless(layer_thickness,icb_heatflux,csb_heatflux,
                     thermal_conductivity)

# Check scalings work - am happy
#if Le <150:
#    with open('scalings/lowLe_Qcmb.pkl', 'rb') as f:
#        slope_Qcmb, intercept_Qcmb = pickle.load(f)
#    with open('scalings/lowLe_den.pkl', 'rb') as f:
#        slope_den, intercept_den = pickle.load(f)
#    with open('scalings/lowLe_prime.pkl', 'rb') as f:
#        [mp1,mp0,ep1,ep0] = pickle.load(f)        
#    with open('scalings/lowLe_v.pkl', 'rb') as f:
#        slope_v, intercept_v = pickle.load(f)        
#else:
#    with open('scalings/highLe_Qcmb.pkl', 'rb') as f:
#        slope_Qcmb, intercept_Qcmb = pickle.load(f)
#    with open('scalings/highLe_den.pkl', 'rb') as f:
#        slope_den, intercept_den = pickle.load(f)
#    with open('scalings/highLe_prime.pkl', 'rb') as f:
#        [mp1,mp0,ep1,ep0] = pickle.load(f)           
#    with open('scalings/highLe_v.pkl', 'rb') as f:
#        slope_v, intercept_v = pickle.load(f) 

#x = Pe*St/Le
#v_est = np.log(slope_v*x+intercept_v) 
##v_est = np.log(m2_v*x**2+m1_v*x+m0_v) 
##v_est=icb_speed
v_est=0
#v_diff = (icb_speed - v_est)/icb_speed*100
#print('v_diff is {:.1f}%'.format(v_diff))
#slope_F = mp1*Pe+mp0
#intercept_F =ep1*Pe+ep0
#F_est = np.log(slope_F*x+intercept_F)
##F_est = np.log(m2_F*x**2+m1_F*x+m0_F) 
F_est = 0
#F_diff = (F-F_est)/F*100
#print('F_diff is {:.1f}%'.format(F_diff))
#densityJump_est = slope_den*x+intercept_den
#Qc_est = slope_Qcmb*x+intercept_Qcmb

# Notes from 25/02/2020
# djdr at CSB
n = profiles.z.size
scale_j = density_solidFe*freezing_speed
Tl = liquidus(csb_radius)
radius = np.linspace(icb_radius,csb_radius,n)
r = np.linspace(icb_radius/csb_radius,1,n)
gravity_prem = premgravity(radius)
density_prem = premdensity(radius)
g0 = gravity_prem[-1]
rho0 = density_prem[-1]
# non-dimensionalise
gravity_prem = gravity_prem/g0
density_prem = density_prem/rho0
grad_gravity_prem = np.gradient(gravity_prem,r)
grad_density_prem = np.gradient(density_prem,r)
# CSB temp
deltaV_liquidFeO_solidFe =getchangevolmelting(8, 8, rho0)
Rrho = rho0/density_solidFe
Rvol = deltaV_solidFe_liquidFe/deltaV_liquidFeO_solidFe
T1 = Tl*heat_capacity*Rrho/(St*latent_heat)


# Variaion of parameters
A = L1*Rrho/(L2*Pe*St*Rvol)
q = Pe*(icb_speed/freezing_speed)/Le + 2/r -2/r # geometry is negligible
q_est = Pe*(v_est/freezing_speed)/Le + 2/r -2/r
c2 = freezing_speed/icb_speed*np.exp(q[-1])
c2_est = freezing_speed/v_est*np.exp(q_est[-1])
c1 = T1 - c2*np.exp(-q[-1])
c1_est = T1 - c2_est*np.exp(-q_est[-1])
T_hom = np.mean(c1+c2*np.exp(-q*r))
T_hom_est = np.mean(c1_est+c2_est*np.exp(-q_est*r))
B=(F*csb_radius/layer_thickness)/T_hom
B_est=(F_est*csb_radius/layer_thickness)/T_hom_est
djdr = A*B*np.exp(F*(csb_radius*r-icb_radius)/layer_thickness) # -2*profiles.solidflux*profiles.oxygen/scale_xi/scale_j
djdr_est = A*B_est*np.exp(F_est*(csb_radius*r-icb_radius)/layer_thickness)
j=cumtrapz(djdr,r,initial=0)
j=0

# Variation of parameters (particular solution)
g = -Pe/(St*Le)*(djdr+2*j/r)
g_est = -Pe/(St*Le)*(djdr_est+2*j/r)
y1 = 1
y2 = np.exp(-q*r)
y2_est = np.exp(-q_est*r)
q_prime = np.gradient(q,r)
q_prime_est = np.gradient(q_est,r)
wronskian = -(q_prime*r+q)*np.exp(-q*r)
wronskian_est = -(q_prime_est*r+q_est)*np.exp(-q_est*r)
u1_prime = -y2*g/wronskian
u2_prime = y1*g/wronskian
u1 = cumtrapz(u1_prime,r,initial=0)
u2 = cumtrapz(u2_prime,r,initial=0)

u1_prime_est = -y2_est*g_est/wronskian_est
u2_prime_est = y1*g_est/wronskian_est
u1_est = cumtrapz(u1_prime_est,r,initial=0)
u2_est = cumtrapz(u2_prime_est,r,initial=0)
# Particular solution
Tp = y1*u1+y2*u2 # particular solution
Tp_prime = np.gradient(Tp,r)
Tp_est = y1*u1_est+y2_est*u2_est # particular solution
Tp_prime_est = np.gradient(Tp_est,r)

# General solution
c2 = freezing_speed/icb_speed*np.exp(q[-1])*(1+Le*Tp_prime[-1]/Pe) # CSB HF
#c2 = freezing_speed*np.exp(q[0]*r[0])*(density_solidFe/(St*density_liquidFe) \
#                           +Le*Tp_prime[0]/Pe)/icb_speed # ICB HF
c1 = T1 - c2*np.exp(-q[-1]) - Tp[-1]
Tc = c1+c2*np.exp(-q*r) # complementary solution 
T = Tc + Tp 

c2_est = freezing_speed/v_est*np.exp(q_est[-1])*(1+Le*Tp_prime_est[-1]/Pe) # CSB HF
#c2 = freezing_speed*np.exp(q[0]*r[0])*(density_solidFe/(St*density_liquidFe) \
#                           +Le*Tp_prime[0]/Pe)/icb_speed # ICB HF
c1_est = T1 - c2_est*np.exp(-q_est[-1]) - Tp_est[-1]
Tc_est = c1_est+c2_est*np.exp(-q_est*r) # complementary solution 
#Tc_est = c1+c2*np.exp(-q_est*r) # complementary solution 
T_est = Tc_est + Tp_est 

# Analytic - Complementary vs full solution
plt.figure()
plt.plot(r,Tc,label = r'$T_c$')
#plt.plot(r,Tp)
plt.plot(r,T,label=r'$T$')
plt.plot(r,Tc_est,label = r'$T_c$ (est.)')
plt.xlabel('r')
plt.legend()

# Plot temperature
plt.figure()
scale_temp = csb_heatflux*1e12/(4*np.pi*csb_radius**2*rho0\
                                *heat_capacity*freezing_speed)
x = pd.to_numeric(profiles.z)/csb_radius
y = pd.to_numeric(profiles.temp)
plt.plot(r,(T)*scale_temp,label='reduced')
plt.plot(x,y,label='numeric')
plt.plot(r,(T_est)*scale_temp,label='estimate')
#plt.plot(x*1e3/csb_radius,y/scale_temp)
plt.xlabel('r')
plt.ylabel(r'$T$')
plt.legend()

# Plot difference between numerical and analytical
difference = ((T)*scale_temp-y)/y*100
print('Maximum percentage difference is {:.2f}%'.format(np.max(np.abs(difference))))
#plt.figure()
#plt.plot(r*csb_radius*1e-3, difference)

# Plot numerical djdr
# fitting
z = np.polyfit(r,np.log(djdr),1)
p = np.poly1d(z)
plt.figure()
plt.plot(r,djdr,label='reduced')
plt.plot(x,profiles.j_grad*csb_radius/scale_j,label='numeric')
#plt.plot(r,np.exp(p(r)),label='fit')
plt.xlabel('r')
plt.ylabel(r'$d\hat{\jmath}/dr$')
#plt.yscale('log')
plt.legend()


# Plot numerical j
j = cumtrapz(djdr,r,initial=0)
j0 = j[-1]
j = j - j0
plt.figure()
plt.plot(r,j,label='reduced')
plt.plot(x,profiles.solidflux/scale_j,label='numeric')
plt.xlabel('r')
plt.ylabel(r'$\hat{\jmath}$')
plt.legend()


# Plot dTdr
dTdr = np.gradient(T,r)
plt.figure()
plt.plot(r,dTdr,label='reduced')
plt.plot(x,profiles.temp_grad*csb_radius/scale_temp,label='numeric')
plt.xlabel('r')
plt.ylabel(r'$d\hat{T}/dr$')
plt.legend()

# Plot dxidr
dxidr = -(L1*gravity_prem*density_prem+dTdr/T)*Rrho/(L2*St*T)
scale_xi,acore = getcsbmassoxygen(pd.to_numeric(inputs.oxygen_bulk)[0],pd.to_numeric(inputs.siliconSulphur_bulk)[0])
plt.figure()
plt.plot(r,dxidr,label='reduced')
plt.plot(x,profiles.xi_grad*csb_radius/scale_xi,label='numeric')
plt.axhline(0,color='k', linestyle='--')
plt.xlabel('r')
plt.ylabel(r'$d\hat{\xi}/dr$')
plt.legend()

# Plot xi
xi = cumtrapz(dxidr,r,initial=0)
xi0 = xi[-1]
xi = xi-xi0+1
plt.figure()
plt.plot(r,xi*scale_xi*acore/aO*100,label='reduced')
plt.plot(x,profiles.oxygen*acore/aO*100,label='numeric')
plt.xlabel('r')
plt.ylabel(r'$\xi$')
plt.legend()

# Plot density
density,_,temp_fluc,xi_fluc,phi_fluc,density_fluc = slurrydensity(profiles.z.to_numpy(),
                                             profiles.temp.to_numpy(),
                                             profiles.oxygen.to_numpy(),
                                             profiles.solidflux.to_numpy(),
                                             layer_thickness,
                                             8, 8, 1e-2) # numeric
density_reduced,_,temp_fluc_reduced,xi_fluc_reduced,phi_fluc_reduced, \
    density_fluc_reduced = slurrydensity(r*csb_radius,T*scale_temp,xi*scale_xi,
                                             j*scale_j,layer_thickness,
                                             8, 8, 1e-2) # reduced
plt.figure()
#plt.plot(r,density_fluc_reduced,label='reduced')
plt.plot(x,density_fluc,label='numeric')
plt.xlabel('r')
plt.ylabel(r'$\rho^\prime$')
plt.legend()

plt.figure()
drho_fluc = np.gradient(density_fluc,x)
dtemp_fluc = np.gradient(temp_fluc,x)
dxi_fluc = np.gradient(xi_fluc,x)
dphi_fluc = np.gradient(phi_fluc,x)
#plt.plot(r,np.gradient(density_fluc_reduced,r),label='reduced')
#colors=plt.cm.GnBu(np.linspace(0.4,1,4))
plt.axhline(0,color='k', linestyle='--')
plt.plot(x*csb_radius*1e-3,dtemp_fluc,color='darkred',lw=2,label=r'$-\rho_{sl} \alpha_T \mathrm{d}T^\prime$')
plt.plot(x*csb_radius*1e-3,dxi_fluc,color='royalblue',lw=2,label=r'$-\rho_{sl} \alpha_\xi \mathrm{d}\xi^\prime$')
plt.plot(x*csb_radius*1e-3,dphi_fluc,color='peru',lw=2,label=r'$\rho_{sl} (\alpha_\phi + \alpha_\xi \xi)\mathrm{d}\phi^\prime$')
plt.plot(x*csb_radius*1e-3,drho_fluc,color='darkgrey',lw=2,label=r'$\mathrm{d} \rho^\prime = \rho_{sl} (-\alpha_T \mathrm{d}T^\prime -  \alpha_\xi \mathrm{d}\xi^\prime +  (\alpha_\phi + \alpha_\xi \xi)\mathrm{d}\phi^\prime)$')
#plt.plot(x,np.gradient(temp_fluc+xi_fluc+phi_fluc,x),label='check',ls='--')
plt.xlabel('Radius (km)')
#plt.ylabel(r'$d \rho^\prime$')
plt.legend(fontsize=11.5)

plt.savefig('figures/eos/'+saveName+'.pdf',format='pdf',dpi=200)

#plt.figure()
#plt.plot(r*csb_radius*1e-3,density_reduced,label='reduced')
#plt.plot(x*csb_radius*1e-3,density,label='numeric')
#radius_prem=np.linspace(icb_radius,csb_radius)
#density_prem=premdensity(radius_prem)
#plt.plot(radius_prem*1e-3,density_prem,'k--')
#plt.xlabel('r')
#plt.ylabel(r'$\rho$')
#plt.legend()