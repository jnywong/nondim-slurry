#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 10:59:46 2018

@author: jennywong
"""

import numpy as np
import pandas as pd

# Profiles
def saveprofiles(outputDir,z,temp,xi,j,phi,density,temp_grad,xi_grad,j_grad, \
                 temp_fluc,xi_fluc,phi_fluc,density_fluc):
    d={'z':z,
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

def saveoutputs(outputDir,F,snow_speed,Q_cmb, Qs, Qs_slurry, 
                     Qs_oc, Ql, Qg, Qg_oc, Qg_slurry, cooling_rate_out, 
                     cmb_temp,acore,state):
    outputs={'F':F,
             'snowSpeed':snow_speed,
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
    
    