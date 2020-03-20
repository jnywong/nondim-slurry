#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 12:02:14 2019

@author: wong
"""

# %% IMPORT

import os
import numpy as np
from scipy.integrate import solve_bvp
#from scipy.special import expi
import pandas as pd

from refparameters import getcsbtemp,geticage, \
getcsbmassoxygen,getsnowspeed,getfreezingspeed,getL1,getL2, \
getLewis,getStefan,getPeclet,getcsbradius,getchangevolmelting
from coreproperties import icb_radius,density_solidFe, \
deltaV_solidFe_liquidFe, heat_capacity,latent_heat
from postprocess import slurrydensity, heatflux
from savescripts import saveprofiles, saveinputs, saveoutputs
from lookup import premgravity, premdensity

def solveslurry(layer_thickness, icb_heatflux, csb_heatflux, thermal_conductivity, \
                ic_age, self_diffusion,sedimentation_constant, \
                mol_conc_oxygen_bulk, mol_conc_SSi, \
                initial_F, n, h):
    
    # %% DIMENSIONLESS NUMBERS
    csb_radius=getcsbradius(layer_thickness)
    mass_conc_O,acore=getcsbmassoxygen(mol_conc_oxygen_bulk,mol_conc_SSi)
    init_snow_speed=getsnowspeed(ic_age) # initial guess
    freezing_speed=getfreezingspeed(icb_heatflux)
    initial_speed = (init_snow_speed + freezing_speed)/freezing_speed # dimensionless
    
    L1,csb_gravity,density0=getL1(csb_radius)
    L2=getL2(mass_conc_O)
    St=getStefan(icb_heatflux,csb_heatflux,csb_radius)
    deltaV_liquidFeO_solidFe=getchangevolmelting(mol_conc_oxygen_bulk,mol_conc_SSi,density0)        
    Le=getLewis(thermal_conductivity,self_diffusion,density0)
    Pe=getPeclet(freezing_speed,csb_radius,self_diffusion)    
    Rrho = density0/density_solidFe
    Rvol = deltaV_solidFe_liquidFe/deltaV_liquidFeO_solidFe
    
    # %% OUTPUT DIRECTORY
    str1=str(np.round(L1,2))
    str2=str(np.round(L2,2))
    str3=str(np.round(St,2))   
    str4=str(np.round(Le,2))
    str5=str(np.round(Pe,2))
    
    if thermal_conductivity==100:
        folderName="highLe/"
    else:
        folderName="lowLe/"
    
    outputDir="results/"+folderName+"L1" + str1 + "_L2" + str2 + "_Le" + str4 + "_St" + str3 + "_Pe" + str5 + "/" 
        
    # Make directory if it doesn't exist
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
#    else:
#        return outputDir
    
    # Ignore cases where csb heat flux is smaller than icb heat flux (actually these are valid solutions)
#    if csb_heatflux < icb_heatflux:
#        return outputDir
    
#    print('(Pe,St) = ({:.2f},{:.2f})'.format(Pe,St))
    
    # Load previous solution to initialise
    if thermal_conductivity==100.0:
        folderName='highLe/'
    else:
        folderName='lowLe/'
    
#    # Previous Stefan number
    m=1
    state=2
    initOn=0
    while state == 2:
        try:
            csb_heatflux_old=csb_heatflux-m*h
            St_old=getStefan(icb_heatflux,csb_heatflux_old,csb_radius)
            str3=str(np.round(St_old,2))        
            inputDir="results/"+folderName+"L1" + str1 + "_L2" + str2 + "_Le" + str4 + "_St" + str3   + "_Pe" + str5 + "/"        
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
                    freezing_speed_old=getfreezingspeed(icb_heatflux_old)
                    Pe_old=getPeclet(freezing_speed_old,csb_radius,self_diffusion)
                    St_old=getStefan(icb_heatflux_old,csb_heatflux,csb_radius)
                    str3=str(np.round(St_old,2)) 
                    str5=str(np.round(Pe_old,2))       
                    inputDir="results/"+folderName+"L1" + str1 + "_L2" + str2 + "_Le" + str4 + "_St" + str3   + "_Pe" + str5 + "/" 
                    data_output=pd.read_csv(inputDir+"outputs.csv",index_col=False)
                    state = np.float(data_output['state'])
                    m=m+1 
                    initOn=2                   
                except FileNotFoundError:
                    initOn=0
                    break
            break
    
#    # Previous L1
#    m=1
#    if layer_thickness == 250e3:
#        state=2
#        initOn=0
#        while state == 2:
#            layer_thickness_old = layer_thickness - 50e3
#            csb_radius_old = icb_radius + layer_thickness_old
#            L1_old,csb_gravity,density0 = getL1(csb_radius_old)
#            freezing_speed_old=getfreezingspeed(icb_heatflux)
#            Pe_old=getPeclet(freezing_speed_old,csb_radius_old,self_diffusion,deltaV_liquidFeO_solidFe)
#            St_old=getStefan(icb_heatflux,csb_heatflux,csb_radius_old,freezing_speed)                    
#            str1=str(np.round(L1_old,2))  
#            str3=str(np.round(St_old,2))   
#            str5=str(np.round(Pe_old,2))
#            inputDir="results/"+folderName+"L1" + str1 + "_L2" + str2 + "_Le" + str4 + "_St" + str3   + "_Pe" + str5 + "/"                
#            data_output=pd.read_csv(inputDir+"outputs.csv",index_col=False)        
#            state = np.float(data_output['state'])
#            initOn=2
    
#    # Read csv data
    if initOn!=0:
        print('Initialised with {}'.format(inputDir))
        data_profile=pd.read_csv(inputDir+"profiles.csv",index_col=False)
        temp0=np.array(data_profile['temp'])
        xi0=np.array(data_profile['oxygen'])
        solid_flux0=np.array(data_profile['solidflux'])
        temp_grad0=np.array(data_profile['temp_grad'])
        initial_F=np.float(data_output['F'])
        snow_speed=np.float(data_output['snowSpeed'])    
        n=temp0.size   

    # %% MESH
    r=np.linspace(icb_radius/csb_radius,1,n) # layer height (m)
    y=np.zeros((4,r.size)) # pre-allocate array

    # %% BOUNDARY VALUE PROBLEM                                         
    def fun_sph(r,y,p): # use ICB or snow speed?
        # Eigenvalues
        F=p[0]
        speed=p[1]
        # PREM gravity and density in layer
        density_prem=premdensity(r*csb_radius)/density0 # dimensionless density
        gravity_prem=premgravity(r*csb_radius)/csb_gravity # dimensionless gravity
        density_grad_prem=np.gradient(density_prem,r)
        gravity_grad_prem=np.gradient(gravity_prem,r)

        term=gravity_prem*density_prem*np.exp(F*(csb_radius*r-icb_radius)/layer_thickness)* \
            (F*csb_radius/layer_thickness+2/r-y[3]/y[0]+ \
             gravity_grad_prem/gravity_prem+density_grad_prem/density_prem)/y[0]
        eq1=-(L1*density_prem*gravity_prem+y[3]/y[0])*Rrho/(L2*St*y[0]) # liquidus
        eq2=(L1*Rrho*term/(L2*St*Pe*Rvol) - eq1*(Rrho*speed+y[2]) - \
             2*y[1]*y[2]/r)/y[1] # oxygen/ dj/dr
        eq3=-Pe/Le*((eq2+2/r*y[2])/St + \
            (speed+2*Le/(r*Pe))*y[3]) # temperature/ d2T/dr2

        return np.vstack([y[3], # eq (15)
                        eq1, # eq (16)
                        eq2, # eq (17)
                        eq3 ]) # eq (18)
                                         
    # Define boundary conditions for solve_bvp
    csb_temp=getcsbtemp(layer_thickness)
#    csb_temp = 6000
  
    def bc_sph(ya,yb,p):
        speed = p[1]

        return np.array([
#                yb[0]-csb_temp*heat_capacity*Rrho/(St*latent_heat),
                yb[0]-csb_temp*heat_capacity*Rrho/(St*latent_heat),
                yb[1]-1, # CSB xi
                ya[2]+speed, # ICB solid flux
                yb[2], # CSB solid flux
                ya[3]+Pe/(St*Le), # ICB heat flux
                yb[3]+Pe/Le # CSB heat flux
                ])

    # Define initial conditions for solve_bvp    
    def ic_sph(y):
        y[0,:]=csb_temp*heat_capacity/(St*latent_heat)*Rrho, # temperature
        y[1,:]=1, # oxygen
        y[2,:]=-initial_speed, # solid flux
        y[3,:]=-Pe/Le # temp gradient
        return y
    
    # Define initial conditions using previous solution
    def ic_old(y):
        y[0,:]=np.reshape(temp0,(1,n)) # temperature
        y[1,:]=np.reshape(xi0,(1,n)) # oxygen
        y[2,:]=np.reshape(solid_flux0,(1,n)) # solid flux
        y[3,:]=np.reshape(temp_grad0,(1,n))
        
        return y
    
    # Non-dimensionalise initial guess
    if initOn!=0:
        temp0=temp0/(csb_heatflux_old*1e12/(4*density0*heat_capacity* \
                    freezing_speed*np.pi*csb_radius**2))
        xi0=xi0/mass_conc_O
        solid_flux0=solid_flux0/(freezing_speed*density_solidFe)
        temp_grad0=temp_grad0/(csb_heatflux_old*1e12/(4*density0*heat_capacity* \
                    freezing_speed*np.pi*csb_radius**3))
#    elif initOn==2:
#        temp0=temp0/csb_temp
#        xi0=xi0/mass_conc_O
#        solid_flux0=solid_flux0/(freezing_speed_old*density0)
#        temp_grad0=temp_grad0/(csb_heatflux*1e12/(4*density0*heat_capacity* \
#                    freezing_speed_old*np.pi*csb_radius**3))
    
    # Run solver       
    tolerance=1e-3
    nmax = 20000
    if initOn!=0:
        sol=solve_bvp(fun_sph,bc_sph,r,ic_old(y),p=[initial_F,initial_speed],tol=tolerance,verbose=2,max_nodes=nmax)
    elif initOn==0:
        print('No initialisation')
        sol=solve_bvp(fun_sph,bc_sph,r,ic_sph(y),p=[initial_F,initial_speed],tol=tolerance,verbose=2,max_nodes=nmax)
    # default tolerance is tol=1e-3
    # default max nodes is 1000 
    
    # If initialisation goes wrong then go for original
    if initOn!=0 and sol.status!=0: # or sol.p[0]<0 or np.log(sol.p[0])>2):
        print('Status = {} - Try again without initialisation'.format(sol.status))
        sol=solve_bvp(fun_sph,bc_sph,r,ic_sph(y),p=[initial_F,snow_speed],tol=tolerance,verbose=2,max_nodes = nmax)
        
    if sol.status==2:
        state=2 # singular Jacobian encountered
        print("Singular Jacobian encountered")
    elif sol.status==1:
        state=1 # max mesh nodes exceeded
        print("Maximum number of mesh nodes exceeded")
    else:
        state=0 # converged
        print("Converged")
        
    # %% OUTPUT
    F_out=sol.p[0]
    print("Barodiffusion parameter is %.2f" % F_out) # F parameter
    icb_speed_out=sol.p[1]*freezing_speed
    snow_speed_out = icb_speed_out - freezing_speed
    
    ic_age_out=geticage(icb_speed_out)
#    print("Snow speed is %.2e m/s" % snow_speed_out) # snow speed
    print("ICB speed is %.2e m/s" % icb_speed_out) # icb speed
    print("IC age is %.2f billion years" % ic_age_out) # ic age
    
    z_out = sol.x*csb_radius
    temp_out = sol.y[0,:]*csb_heatflux*1e12/(4*density0*heat_capacity* \
                    freezing_speed*np.pi*csb_radius**2)
    xi_out= sol.y[1,:]*mass_conc_O
    j_out= sol.y[2,:]*freezing_speed*density_solidFe
    j_out[-1]=0. # BVP minimises j to be close to, but not exactly zero at the CSB
    temp_grad_out= sol.yp[0,:]*csb_heatflux*1e12/(4*density0*heat_capacity* \
                    freezing_speed*np.pi*csb_radius**3)
    xi_grad_out= sol.yp[1,:]*mass_conc_O/csb_radius
    j_grad_out=sol.yp[2,:]*freezing_speed*density_solidFe/csb_radius
    
    # %% POST-PROCESS
    # Slurry density
    density_slurry,phi_out,temp_fluc,xi_fluc,phi_fluc,density_fluc=slurrydensity(z_out,temp_out,xi_out,j_out, \
                                         layer_thickness,mol_conc_oxygen_bulk, \
                                         mol_conc_SSi,sedimentation_constant)
    density_jump=density_slurry[0]-density_slurry[-1]
    print("Density jump is {:.2f} kg/m^-3".format(density_jump))
    
    # Heat balance
    (Q_cmb, Qs, Qs_slurry, Qs_oc, Ql, Qg, Qg_oc, Qg_slurry, cooling_rate_out,
            cmb_temp, temp_ad) = heatflux(z_out,temp_out,xi_out,j_out,phi_out, \
                             temp_grad_out,xi_grad_out, density_slurry, \
                             snow_speed_out,freezing_speed,icb_heatflux*1e12,
                             layer_thickness,thermal_conductivity,csb_heatflux*1e12,n)

    # %% SAVE 
    saveprofiles(outputDir,z_out,temp_out,xi_out,j_out,phi_out,density_slurry, \
                 temp_grad_out,xi_grad_out,j_grad_out,temp_fluc,xi_fluc, \
                 phi_fluc,density_fluc)
    saveinputs(outputDir,n,layer_thickness,thermal_conductivity,icb_heatflux, \
           csb_heatflux, mol_conc_oxygen_bulk,mol_conc_SSi, \
           self_diffusion,sedimentation_constant)
    saveoutputs(outputDir,F_out,snow_speed_out,Q_cmb, Qs, Qs_slurry, \
                     Qs_oc, Ql, Qg, Qg_oc, Qg_slurry, cooling_rate_out, \
                     cmb_temp,acore,state)
    print('Run ' + outputDir[0:-1] + ' is saved')
    
    return (outputDir, z_out, temp_out, xi_out, j_out, density_slurry, csb_temp)