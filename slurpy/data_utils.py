#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:10:46 2020

@author: jennywong
"""
import numpy as np
import pandas as pd

import slurpy.getparameters as gp
from slurpy.coreproperties import density_solidFe, deltaV_solidFe_liquidFe

# Profiles
def saveprofiles(outputDir,z,temp,xi,j,phi,density,temp_grad,xi_grad,j_grad, \
                 temp_fluc,xi_fluc,phi_fluc,density_fluc):
    d={'r':z,
       'temp':temp,
       'oxygen':xi,
       'solidflux':j,
       'phi':phi,
       'density':density,
       'temp_grad':temp_grad,
       'xi_grad':xi_grad,
       'j_grad':j_grad,
       'temp_fluc':temp_fluc,
       'xi_fluc':xi_fluc,
       'phi_fluc':phi_fluc,
       'density_fluc':density_fluc}
    df=pd.DataFrame(d)
    df.to_csv(outputDir+'profiles.csv',index=False)

def saveinputs(outputDir,n,layer_thickness,thermal_conductivity,icb_heatflux, \
                    csb_heatflux,mol_conc_oxygen_bulk,mol_conc_SSi, \
                    self_diffusion,sedimentation_constant):
    inputs={'n':n,'layerThickness':layer_thickness,
            'thermalConductivity':thermal_conductivity,
            'icb_heatflux':icb_heatflux,
            'csb_heatflux':csb_heatflux,
            'oxygen_bulk':mol_conc_oxygen_bulk,
            'siliconSulphur_bulk':mol_conc_SSi,
            'selfdiffusionCoefficient':self_diffusion,
            'sedimentationConstant':sedimentation_constant
            }
    df=pd.DataFrame(inputs,index=[0]) # scalar data requires passing index=[0]
    df.to_csv(outputDir+'inputs.csv')

def saveoutputs(outputDir,F,snow_speed,icb_speed, Q_cmb, Qs, Qs_slurry, 
                     Qs_oc, Ql, Qg, Qg_oc, Qg_slurry, cooling_rate_out, 
                     cmb_temp,acore,state):
    outputs={'F':F,
             'snowSpeed':snow_speed,
             'icbSpeed':icb_speed,             
             'Q_cmb':Q_cmb,
             'Qs':Qs,
             'Qs_slurry': Qs_slurry,
             'Qs_oc': Qs_oc,
             'Ql': Ql,
             'Qg': Qg,
             'Qg_oc':Qg_oc,
             'Qg_slurry': Qg_slurry,
             'cooling_rate_out': cooling_rate_out,
             'cmb_temp': cmb_temp,
             'a_core':acore,
             'state':state}
    df=pd.DataFrame(outputs,index=[0]) # scalar data requires passing index=[0]
    df.to_csv(outputDir+'outputs.csv')
    
def savephaseboundaries(outputDir,nPe,nSt,Pe,St,density_jump,cmb_heatflux,
                            delta, F, icb_speed):
    nPe=np.array(nPe,ndmin=1)
    np.savetxt(outputDir+'nPe.csv',nPe,fmt='%i')
    nSt=np.array(nSt,ndmin=1)
    np.savetxt(outputDir+'nSt.csv',nSt,fmt='%i')
    np.savetxt(outputDir+'Pe.csv',Pe)
    np.savetxt(outputDir+'St.csv',St)
    np.savetxt(outputDir+'density_jump.csv',density_jump)
    np.savetxt(outputDir+'cmb_heatflux.csv',cmb_heatflux)
    np.savetxt(outputDir+'delta.csv',delta)
    np.savetxt(outputDir+'F.csv',F)
    np.savetxt(outputDir+'icb_speed.csv',icb_speed)    

def readdata(inputDir):
    # Data from spherical model
    data_inputs=pd.read_csv(inputDir+"/inputs.csv",index_col=False)
    data_outputs=pd.read_csv(inputDir+"/outputs.csv",index_col=False)
    data_profiles=pd.read_csv(inputDir+"/profiles.csv",index_col=False)
    
    return data_inputs,data_outputs,data_profiles    

def get_outputDir(layer_thickness,icb_heatflux,csb_heatflux,thermal_conductivity,
                  mol_conc_oxygen_bulk=8.,self_diffusion=0.98e-8,mol_conc_SSi=8,
                  model='prem'):
    csb_radius=gp.getcsbradius(layer_thickness)
    mass_conc_O,acore=gp.getcsbmassoxygen(mol_conc_oxygen_bulk)
    freezing_speed=gp.getfreezingspeed(icb_heatflux)
    Lip,csb_gravity,density0=gp.getLip(csb_radius,model)
    Lix=gp.getLix(mass_conc_O)
    St=gp.getStefan(icb_heatflux,csb_heatflux,csb_radius)
    Le=gp.getLewis(thermal_conductivity,self_diffusion,density0)
    Pe=gp.getPeclet(freezing_speed,csb_radius,self_diffusion)
    
    str1=str(np.round(Le,2)).replace('.','_')
    str2=str(np.round(Lip,2)).replace('.','_')
    str3=str(np.round(Lix,2)).replace('.','_')
    str4=str(np.round(Pe,2)).replace('.','_')
    str5=str(np.round(St,2)).replace('.','_')

    foldername = "Le_{}".format(str1)
    filename = "Lip_{}_Lix_{}_Pe_{}_St_{}".format(str2,str3,str4,str5)
    
    return foldername,filename