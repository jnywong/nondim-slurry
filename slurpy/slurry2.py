#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 13:20:26 2020

@author: Jenny Wong

###############################################################################
# SLURRY.PY                                                                   #
###############################################################################

# Solve the 1D, steady, spherical slurry system outlined in Wong et al. (in prep)
# (see also Wong et al. 2018).

Parameters
----------
plotOn : int
    Toggle 0='off' and 1='on' to display the solution output. Recommended off for large parameter searches but useful for probing a small number of solutions
layer_thicknesses : numpy array
    Specify slurry layer thickness (m).
thermal_conductivities : numpy array
    Specify thermal conductivity (W/m/K).
icb_heatfluxes : numpy array
    Specify heat flux through the ICB (TW).
csb_heatfluxes : numpy array
    Specify heat flux through the CSB (TW).
h : float
    Specify step size through the heat fluxes for parameter searches.
sensitivityOn : int
    Toggle manual control of the CSB temperature and/or CSB oxygen concentration for sensitivity studies.
mol_conc_oxygen_bulk : float
    If sensitivityOn=1, then manually specify the CSB oxygen concentration. Default value is 8 mol.%.
csb_temp : float
    If sensitivityOn=1, then manually specify the CSB temperture (K). Default value depends on the layer thickness, given by function slurpy.lookup.liquidus.
sedimentation_constant : float
    If sensitivityOn=1, then manually specify the sedimentation constant used to derive the solid fraction from the solid flux. Default value is 1e-2 kg s/m^3.

"""

import numpy as np
import os
import pandas as pd
import slurpy.getparameters as gp
# import slurpy.coreproperties as cp
import slurpy.lookup as lp
import slurpy.data_utils as sp

from scipy.integrate import solve_bvp
from slurpy.coreproperties import icb_radius, deltaV_solidFe_liquidFe, \
    density_solidFe, heat_capacity, latent_heat, year

def solveslurry(layer_thickness, icb_heatflux, csb_heatflux, thermal_conductivity, \
            csb_temp, h, mol_conc_oxygen_bulk=8, sedimentation_constant=1e-2,
            self_diffusion=0.98e-8, mol_conc_SSi=8, model = 'prem', \
            initial_F=5, initial_icAge=0.5, maxSt=6, n=100,
            tolerance=1e-3,nmax=2e4):

    # %% Nested functions to define BVP
    def fun(r,y,p):
        # Eigenvalues
        F=p[0]
        speed=p[1]
        # PREM gravity and density in layer
        if model =='prem':
            density_seis=lp.premdensity(r*csb_radius)/density0 # dimensionless density
        elif model == 'ohtaki':
            _,density_seis = lp.ohtaki(r*csb_radius)/density0
        gravity_seis=lp.premgravity(r*csb_radius)/csb_gravity # dimensionless gravity
        density_grad_seis=np.gradient(density_seis,r)
        gravity_grad_seis=np.gradient(gravity_seis,r)

        # r derivative of barodiffusion term
        term=gravity_seis*density_seis*np.exp(F*(csb_radius*r-icb_radius)/layer_thickness)* \
            (F*csb_radius/layer_thickness+2/r-y[3]/y[0]+ \
             gravity_grad_seis/gravity_seis+density_grad_seis/density_seis)/y[0]
        # liquidus (=dxi/dr)
        eq1=-(Lip*density_seis*gravity_seis+y[3]/y[0])*Rrho/(Lix*St*y[0])
        # oxygen eqn (=dj/dr)
        eq2=(Lip*Rrho*term/(Lix*St*Pe*Rvol) - eq1*(Rrho*speed*(icb_radius/csb_radius)**2+y[2]) - \
             2*y[1]*y[2]/r)/y[1]
        # temp eqn (=d2T/dr2)
        eq3=-Pe/Le*((eq2+2/r*y[2])/St + \
            (speed*(icb_radius/csb_radius)**2+2*Le/(r*Pe))*y[3])
        return np.vstack([y[3],eq1,eq2,eq3])

    def bcs(ya,yb,p):
        speed = p[1]
        return np.array([
                yb[0]-csb_temp*heat_capacity*Rrho/(St*latent_heat), # CSB temp
                yb[1]-1, # CSB xi
                ya[2]+speed, # ICB solid flux
                yb[2], # CSB solid flux
                ya[3]+Pe/(St*Le), # ICB heat flux
                yb[3]+Pe/Le # CSB heat flux
                ])

    def ics(y):
        y[0,:]=csb_temp*heat_capacity*Rrho/(St*latent_heat), # temperature
        y[1,:]=1, # oxygen
        y[2,:]=-initial_speed, # solid flux
        y[3,:]=-Pe/Le # temp gradient
        return y

    # Define initial conditions using previous solution
    def ic_old(y):
        y[0,:]=np.reshape(temp0,(1,n)) # temperature
        y[1,:]=np.reshape(xi0,(1,n)) # oxygen
        y[2,:]=np.reshape(solid_flux0,(1,n)) # solid flux
        y[3,:]=np.reshape(temp_grad0,(1,n)) # temp gradient
        return y

    # %% DIMENSIONLESS NUMBERS
    csb_radius=gp.getcsbradius(layer_thickness)
    mass_conc_O,acore=gp.getcsbmassoxygen(mol_conc_oxygen_bulk)
    init_snow_speed=gp.getsnowspeed(initial_icAge) # initial guess
    freezing_speed=gp.getfreezingspeed(icb_heatflux)
    initial_speed = (init_snow_speed + freezing_speed)/freezing_speed # dimensionless
    Lip,csb_gravity,density0=gp.getLip(csb_radius,model)
    Lix=gp.getLix(mass_conc_O)
    St=gp.getStefan(icb_heatflux,csb_heatflux,csb_radius)
    deltaV_liquidFeO_solidFe=gp.getchangevolmelting(mol_conc_oxygen_bulk,density0)
    Le=gp.getLewis(thermal_conductivity,self_diffusion,density0)
    Pe=gp.getPeclet(freezing_speed,csb_radius,self_diffusion)
    Rrho=density0/density_solidFe
    Rvol = deltaV_solidFe_liquidFe/deltaV_liquidFeO_solidFe
    scale_temp = gp.get_tempScale(csb_heatflux,csb_radius,density0,freezing_speed)
    scale_xi = mass_conc_O
    scale_j = gp.get_jScale(freezing_speed)


    # %% OUTPUT DIRECTORY
    str1=str(np.round(Le,2)).replace('.','_')
    str2=str(np.round(Lip,2)).replace('.','_')
    str3=str(np.round(Lix,2)).replace('.','_')
    str4=str(np.round(Pe,2)).replace('.','_')
    str5=str(np.round(St,2)).replace('.','_')

    if model == 'ohtaki':
        outputDir="ohtaki/Le_{}/Lip_{}_Lix_{}_Pe_{}_St_{}/".format(str1,str2,str3,str4,str5)
    else:
        outputDir="results/Le_{}/Lip_{}_Lix_{}_Pe_{}_St_{}/".format(str1,str2,str3,str4,str5)
    # Make directory if it doesn't exist
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    # Ignore cases where csb heat flux is smaller than icb heat flux
    elif csb_heatflux < icb_heatflux:
        return (outputDir,0,0,0,0,0)
    # Impose upper limit on St
    elif St> maxSt:
        return (outputDir,0,0,0,0,0)
    # Skip if directory already exists
#    else:
#        return (outputDir,0,0,0,0,0)

    # Load previous solution to initialise
    # Previous Stefan number
    m=1
    state=2
    initOn=0
    # FIX: loop back through all of St before Pe
    while state == 2:
        try:
            csb_heatflux_old=csb_heatflux-m*h
            St_old=gp.getStefan(icb_heatflux,csb_heatflux_old,csb_radius)
            str5=str(np.round(St_old,2)).replace('.','_')
            inputDir="results/Le_{}/Lip_{}_Lix_{}_Pe_{}_St_{}/".format(str1,str2,str3,str4,str5)
            data_output=pd.read_csv(inputDir+"outputs.csv",index_col=False)
            state = np.float(data_output['state'])
            m=m+1
            initOn=1
        except FileNotFoundError:
            # Previous Peclet number
            m=1
            state=2
            while state == 2:
                try:
                    icb_heatflux_old = icb_heatflux-m*h
                    freezing_speed_old=gp.getfreezingspeed(icb_heatflux_old)
                    Pe_old=gp.getPeclet(freezing_speed_old,csb_radius,self_diffusion)
                    St_old=gp.getStefan(icb_heatflux_old,csb_heatflux,csb_radius)
                    str4=str(np.round(Pe_old,2)).replace('.','_')
                    str5=str(np.round(St_old,2)).replace('.','_')
                    inputDir="results/Le_{}/Lip_{}_Lix_{}_Pe_{}_St_{}/".format(str1,str2,str3,str4,str5)
                    data_output=pd.read_csv(inputDir+"outputs.csv",index_col=False)
                    state = np.float(data_output['state'])
                    m=m+1
                    initOn=2
                except FileNotFoundError:
                    initOn=0
                    break
            break

    # Read csv data
    if initOn!=0:
        data_profile=pd.read_csv(inputDir+"profiles.csv",index_col=False)
        temp0=np.array(data_profile['temp'])
        xi0=np.array(data_profile['oxygen'])
        solid_flux0=np.array(data_profile['solidflux'])
        temp_grad0=np.array(data_profile['temp_grad'])
        initial_F=np.float(data_output['F'])
        snow_speed=np.float(data_output['snowSpeed'])
        initial_speed=np.float(data_output['icbSpeed'])
        freezing_speed_old=initial_speed-snow_speed
        n=temp0.size
        # Non-dimensionalise initial guess
        scale_tempOld = gp.get_tempScale(csb_heatflux_old,csb_radius,density0,
                                      freezing_speed)
        scale_jOld=gp.get_jScale(freezing_speed_old)
        temp0=temp0/scale_tempOld
        xi0=xi0/scale_xi
        solid_flux0=solid_flux0/scale_jOld
        temp_grad0=temp_grad0/scale_tempOld*csb_radius


    # %% MESH
    r=np.linspace(icb_radius/csb_radius,1,n)
    y=np.zeros((4,r.size)) # pre-allocate soln array

    # %% BOUNDARY VALUE PROBLEM
    # Run solver - default solver tolerance is 1e-3, default max nodes is 1000
    if initOn!=0:
        print('Initialised with {}'.format(inputDir))
        sol=solve_bvp(fun,bcs,r,ic_old(y),p=[initial_F,initial_speed],tol=tolerance,verbose=2,max_nodes=nmax)
    elif initOn==0:
        print('No initialisation')
        sol=solve_bvp(fun,bcs,r,ics(y),p=[initial_F,initial_speed],tol=tolerance,verbose=2,max_nodes=nmax)

    # If initialisation gives no soln then try no initialisation
    if initOn!=0 and sol.status!=0:
        print('Status = {} - Try again without initialisation'.format(sol.status))
        sol=solve_bvp(fun,bcs,r,ics(y),p=[initial_F,initial_speed],tol=tolerance,verbose=2,max_nodes = nmax)

    if sol.status==2:
        state=2
        print("Singular Jacobian encountered")
    elif sol.status==1:
        state=1
        print("Maximum number of mesh nodes exceeded")
    else:
        state=0
        print("Converged")

    # %% OUTPUT
    F_out=sol.p[0]
    icb_speed_out=sol.p[1]*freezing_speed
    snow_speed_out = icb_speed_out - freezing_speed
    ic_age_out=gp.geticage(icb_speed_out)
    print("Mixing parameter is %.2f" % F_out)
    print("Freezing speed is {:.2e} m/s = {:.3f} mm/yr".format(freezing_speed,freezing_speed*1e3*year))
    print("ICB speed is {:.2e} m/s = {:.3f} mm/yr".format(icb_speed_out,icb_speed_out*1e3*year))
    print("IC age is %.2f Ga" % ic_age_out)

    # Nondimensional to dimensional
    r_out = sol.x*csb_radius
    temp_out = sol.y[0,:]*scale_temp
    xi_out= sol.y[1,:]*scale_xi
    j_out= sol.y[2,:]*scale_j
    j_out[-1]=0. # BVP minimises j to be close to, but not exactly zero at the CSB
    temp_grad_out= sol.yp[0,:]*scale_temp/csb_radius
    xi_grad_out= sol.yp[1,:]*scale_xi/csb_radius
    j_grad_out=sol.yp[2,:]*scale_j/csb_radius

    # %% POST-PROCESS
    # Slurry density
    density,phi_out,temp_denFluc,xi_denFluc,phi_denFluc,density_fluc= \
        gp.slurrydensity(r_out,temp_out,xi_out,j_out, \
                         mol_conc_oxygen_bulk,sedimentation_constant,model=model)
    density_jump=density[0]-density[-1]
    print("Density jump is {:.2f} kg/m^-3".format(density_jump))
    # Stable layer?
    density_fluc_grad = np.gradient(density_fluc,r_out)
    unstable = density_fluc_grad[density_fluc_grad>0]
    if unstable.size!=0:
        print('Unstable slurry')
    else:
        print('Stable slurry')

    # Heat balance
    (Q_cmb, Qs, Qs_slurry, Qs_oc, Ql, Qg, Qg_oc, Qg_slurry, cooling_rate_out,
            cmb_temp, temp_ad) = gp.heatflux(r_out,temp_out,xi_out,j_out,phi_out, \
                             temp_grad_out,xi_grad_out, density, \
                             icb_speed_out,icb_heatflux*1e12,
                             layer_thickness,thermal_conductivity,csb_heatflux*1e12,n)

    # %% SAVE
    sp.saveprofiles(outputDir,r_out,temp_out,xi_out,j_out,phi_out,density, \
                 temp_grad_out,xi_grad_out,j_grad_out,temp_denFluc,xi_denFluc, \
                 phi_denFluc,density_fluc)
    sp.saveinputs(outputDir,n,layer_thickness,thermal_conductivity,icb_heatflux, \
           csb_heatflux, mol_conc_oxygen_bulk,mol_conc_SSi, \
           self_diffusion,sedimentation_constant)
    sp.saveoutputs(outputDir,F_out,snow_speed_out,icb_speed_out,Q_cmb, Qs, Qs_slurry, \
                     Qs_oc, Ql, Qg, Qg_oc, Qg_slurry, cooling_rate_out, \
                     cmb_temp,acore,state)
    print('Run {} is saved'.format(outputDir))

    return (outputDir, r_out, temp_out, xi_out, j_out, F_out, icb_speed_out, density)
