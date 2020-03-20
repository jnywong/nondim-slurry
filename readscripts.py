#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:41:26 2019

@author: wong
"""

import pandas as pd

def read_box1data(layer_thickness,icb_heatflux,thermal_conductivity,ic_age):
    str1=str(layer_thickness*1e-3).replace(".0","")
    str2=str(icb_heatflux).replace(".0","")
    str3=str(thermal_conductivity).replace(".0","")
    str4=str(ic_age)
    inputDir="box1_d"+str1+"_icbhf"+str2+"_k"+str3+"_icage"+str4+"/"
    
    data_inputs=pd.read_csv(inputDir+"inputs.csv",index_col=False)
    data_outputs=pd.read_csv(inputDir+"outputs.csv",index_col=False)
    data_profiles=pd.read_csv(inputDir+"profiles.csv",index_col=False)
    
    return data_inputs,data_outputs,data_profiles

def read_box2data(layer_thickness,icb_heatflux,csb_heatflux,thermal_conductivity):
    str1=str(layer_thickness*1e-3).replace(".0","")
    str2=str(icb_heatflux).replace(".0","")
    str3=str(thermal_conductivity).replace(".0","")
#    str4=str(round(csb_heatflux[0],2))
    str4=str(csb_heatflux)
    inputDir="box2_d"+str1+"_icbhf"+str2+"_k"+str3+"_csbhf"+str4+"/"
    
    # Data from box model 2 eigenvalues
    data_inputs=pd.read_csv(inputDir+"inputs.csv",index_col=False)
    data_outputs=pd.read_csv(inputDir+"outputs.csv",index_col=False)
    data_profiles=pd.read_csv(inputDir+"profiles.csv",index_col=False)
    
    return data_inputs,data_outputs,data_profiles

def read_sphdata(inputDir):
    # Data from spherical model
    data_inputs=pd.read_csv(inputDir+"inputs.csv",index_col=False)
    data_outputs=pd.read_csv(inputDir+"outputs.csv",index_col=False)
    data_profiles=pd.read_csv(inputDir+"profiles.csv",index_col=False)
    
#    data_outputs=0
    
    return data_inputs,data_outputs,data_profiles