#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 10:12:15 2019
@author: jennywong
"""

# Scripts to deal with results from parameter search

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lookup import premdensity
from coreproperties import icb_radius, aO
from refparameters import getL1, getL2, getStefan, getLewis, getPeclet, \
    getcsbmassoxygen, getcsbradius, getchangevolmelting, getfreezingspeed

def makeArray(text):
    return np.fromstring(text,sep=' ')

# %% Does this run contain a stable layer? If so, what is the density jump and 
#    CMB heat flux?
def stablelayer(icb_heatflux,layer_thickness,csb_heatflux,thermal_conductivity,plotOn):
    if thermal_conductivity==100.0:
        folderName='highLe/'
    else:
        folderName='lowLe/'

    mol_conc_oxygen_bulk = 8.
    mol_conc_SSi = 8.
    self_diffusion= 0.98e-8
    
    csb_radius=getcsbradius(layer_thickness)
    mass_conc_O, acore = getcsbmassoxygen(mol_conc_oxygen_bulk, mol_conc_SSi)
    density0 = premdensity(csb_radius)
#    deltaV_liquidFeO_solidFe=getchangevolmelting(mol_conc_oxygen_bulk, mol_conc_SSi,density0)
    freezing_speed=getfreezingspeed(icb_heatflux)
    
    L1,_,_=getL1(csb_radius)
    L1=np.float(L1)
    L2=np.float(getL2(mass_conc_O))
    Le=np.float(getLewis(thermal_conductivity,self_diffusion,density0))
    Pe=np.float(getPeclet(freezing_speed,csb_radius,self_diffusion))
    St=np.float(getStefan(icb_heatflux,csb_heatflux,csb_radius))
    
    if icb_heatflux >= csb_heatflux:
        density_jump=0
        cmb_heatflux=0
        stable_density=0
        stable_oxygen=0
        solution=0
        density_exceeds_prem=0
        F=0
        icb_speed=0
        return (density_jump,cmb_heatflux,stable_density,stable_oxygen,solution,
                density_exceeds_prem,L1,L2,St,Pe,Le,F,icb_speed)

    str1=str(np.round(L1,2))
    str2=str(np.round(L2,2))
    str3=str(np.round(St,2))
    str4=str(np.round(Le,2))
    str5=str(np.round(Pe,2))
    outputDir="results/"+folderName+"L1" + str1 + "_L2" + str2 + "_Le" + str4 + "_St" + str3   + "_Pe" + str5 + "/"
    # Read csv data
    data_profile=pd.read_csv(outputDir+"profiles.csv",index_col=False)
    z=np.array(data_profile['z'])
    density=np.array(data_profile['density'])
    xi=np.array(data_profile['oxygen'])
    
    # Check for no solutions
    data_output=pd.read_csv(outputDir+"outputs.csv",index_col=False)
    F=np.array(data_output['F'])
    snow_speed =np.array(data_output['snowSpeed'])
    icb_speed = freezing_speed + snow_speed
    
    state=np.array(data_output['state'])
    if np.isnan(density).any()==True or state!=0:
        print(outputDir+':No solution')
        density_jump=0
        cmb_heatflux=0
        stable_density=0
        stable_oxygen=0
        solution=0
        density_exceeds_prem=0
        F=0
        icb_speed=0
        return (density_jump,cmb_heatflux,stable_density,stable_oxygen,solution,
                density_exceeds_prem,L1,L2,St,Pe,Le,F,icb_speed)
    else:
        solution=1
    
    # Check density fluc and oxygen fluc gradients
    # if drho'/dr < 0 then stable, if d xi'/dr >0 then stable
    radius=z*1e-3
    density_fluc=np.array(data_profile['density_fluc'])
    xi_fluc=np.array(data_profile['xi_fluc'])
    densityGrad=np.gradient(density_fluc,radius)
    a=np.where(densityGrad>0)
    xi_gradient=np.gradient(xi_fluc,radius)
    b=np.where(xi_gradient<0)
    if len(a[0])==0:
        print(outputDir+':stable density gradient')
        stable_density=1
    elif len(a[0])!=0:
        print(outputDir+':unstable density gradient')
        stable_density=0
    if len(b[0])!=0:
        print(outputDir+':stable oxygen gradient')
        stable_oxygen=1
    elif len(b[0])==0: 
        print(outputDir+':unstable oxygen gradient')
        stable_oxygen=0

    # Density jump
    density_jump=density[0]-density[-1]
    # CMB heat flux
    cmb_heatflux= np.array(data_output['Q_cmb'])
                    
    # Does slurry density exceed PREM?
    density_prem=premdensity(radius*1e3)
    c=np.where(density<density_prem)
    if (len(c[0])>1):
        print(outputDir+":slurry density smaller than PREM")
        density_exceeds_prem=0    
    else:
        print(outputDir+":slurry density exceeds PREM")    
        density_exceeds_prem=1
        
    # Is slurry partially stratified?        
    threshold = icb_radius+100e3 # assume seismic data resolution is at least 100km # note z is non-uniform
    idx_threshold = np.argwhere(z<threshold)
    len_threshold = idx_threshold.size
    if (len(a[0])==0):
        print(outputDir+":fully stratified")    
        stable_density=1
    elif (len(a[0])<len_threshold) and density_exceeds_prem==1:
        # if more than 2/3 (thinnest seismic layer/150km) then partially stratified
        print(outputDir+":partially stratified")  
        stable_density=2 # partially stratified
    else:
        print(outputDir+":slurry density gradient smaller than PREM gradient")
        stable_density=0
    
    # Plot density
    if plotOn==1:
        fig1 = plt.figure
        fig1, ax1 = plt.subplots()
        ax1.plot(radius,density)
        ax1.plot(radius,density_prem,'k--')
        ax1.set(xlabel='radius ($km$)',ylabel='density ($kgm^{-3}$)')
        ax1.legend(["slurry","PREM"])
        fig2 = plt.figure
        fig2, ax2 = plt.subplots()
        ax2.plot(radius,densityGrad)
        ax2.hlines(0,radius[0],radius[-1],'k',linestyles= '--')
#        ax2.plot(radius,densityGrad_prem,'k--')
#        ax2.legend(["slurry","PREM"])
        ax2.set(xlabel='radius ($km$)',ylabel='density fluctuation gradient ($kgm^{-4}$)')
        fig3 = plt.figure
        fig3, ax3 = plt.subplots()
        acore=np.array(data_output['a_core'])
        ax3.plot(radius,xi_fluc*acore/aO*100)
        ax3.set(xlabel='radius ($km$)',ylabel='oxygen concentration fluctuations ($mol.\%$)')
        fig4 = plt.figure
        fig4, ax4 = plt.subplots()
        ax4.plot(radius,xi_gradient*acore/aO*100)
        ax4.hlines(0,radius[0],radius[-1],'k',linestyles= '--')
        ax4.set(xlabel='radius ($km$)',ylabel='oxygen concentration gradient ($mol.\% m^{-1}$)')
        plt.show()
    
    return (density_jump,cmb_heatflux,stable_density,stable_oxygen,solution,
            density_exceeds_prem,L1,L2,St,Pe,Le,F,icb_speed)

#%%----------------------------------------------------------------------------   
def fitphaseboundary(Pe,St,a,b,degree,n=100):
    # Polynomial fit
#    line_phase=np.polyfit(Pe,St,degree,full=False) #,full=True)
#    poly_fit=np.poly1d(line_phase)
    xp = np.linspace(a, b, n)
#    yp = poly_fit(xp)
#    yp=spline(Pe,St,xp)
#    tck=interpolate.splrep(Pe,St) #,degree)
#    yp=interpolate.splev(xp,tck)
    z=np.polyfit(Pe,St,degree)
    p=np.poly1d(z)
    yp=p(xp)
    return (xp,yp) #,line_phase)

#%% 
def getphaseboundary(Pe,St,delta,nPe,nSt,nGroup=100):
    # determine the edges of the phase boundaries
    # delta is the matrix containing regime classification for each Pe,St:
    #   0 - density is less than PREM or density fluc gradient >0 -> REJECT
    #   1 - stable density and stable oxygen -> ACCEPT
    #   2 - stable density and unstable oxygen -> ACCEPT
    #   3 - partially stratified and stable oxygen -> ACCEPT
    #   4 - partially startified and unstable oxygen -> ACCEPT
    #   5 - no solution -> REJECT
    
    # delta needs correcting because of numerical artefacts
    delta_correct=np.zeros((nPe,nSt))
    
    # Initialise lists
    # Lower boundary between stable/unstable and no solution
    Pe_lower=[]
    St_lower=[]
    # Upper boundary between stable and unstable/no solution
    Pe_upper=[]
    St_upper=[]
    # Boundary between stable and unstable
    Pe_unstable=[] 
    St_unstable=[]
    # Boundary between stable O and unstable O
    Pe_stable=[] 
    St_stable=[]
    # Vertical boundary between no solution and unstable
    Pe_unstable_vertical=[]
    St_unstable_vertical=[]
    # Vertical boundary between unstable and stable (low Le)
    Pe_stable_vertical=[]
    St_stable_vertical=[]
    
    seq=[]
    for i in range(nPe): 
        m=0
        groupStart=[]
        groupEnd=[]
        groupID=[] # to be filled with classification 0-5
        
        groupStart.append(int(0))
        groupID.append(delta[i,groupStart[m]])
        # Find boundaries
        for j in range(nSt):
            if delta[i,j]!=groupID[m]:
                groupStart.append(j)
                groupEnd.append(j-1)
                groupID.append(delta[i,j])
                m=m+1                
        groupEnd.append(nSt)
        # Remove artifacts from groupID sequence
        sequence_correct=np.asarray(groupID)
        nSeq=sequence_correct.size
        for l in range(nSeq):
            if l!=0 and l!=nSeq-1 and sequence_correct[l-1]==groupID[l+1] and sequence_correct[l-1]!=5 and groupID[l+1]!=5 and groupID[l]==0:
                sequence_correct[l]=sequence_correct[l-1]   
#            elif l==nSeq-2 and groupID[l+1]==5 and groupID[l]==0:
#                sequence_correct[l]=5
            elif l!=0 and l!=nSeq-1 and groupID[l]==5 and groupID[l+1]==1:
                sequence_correct[l]=1
            elif l!=0 and l!=nSeq-1 and sequence_correct[l-1]==5 and groupID[l+1]==5:
                sequence_correct[l]=groupID[l] 
            elif l!=0 and l!=nSeq-1 and groupID[l-1]==1 and groupID[l]==5 and groupID[l+1]==2:
                sequence_correct[l]=1
            elif l!=0 and l!=nSeq-1 and groupID[l-1]!=5 and groupID[l]==5 and groupID[l+1]!=5 and groupID[l+1]!=0:
                sequence_correct[l]=groupID[l-1]
            else:
                sequence_correct[l]=groupID[l]
        # Apply new sequence to delta
        for k in range(nSeq):
            delta_correct[i,groupStart[k]:groupEnd[k]+1]=sequence_correct[k]
            
        delta_correct = delta
        
            
        # Make delta_correct for unstable solutions the same value  
#        delta_correct[i,:] = np.where(np.logical_or(delta_correct[i,:]==5,delta_correct[i,:]==1), delta_correct[i,:], 0)
#        sequence_correct = np.where(np.logical_or(sequence_correct==5,sequence_correct==1), sequence_correct, 0)
#        seq.append(sequence_correct)
        
        # Group phase boundaries #CHANGE TO LOOK AT SEQUENCE AWAY FROM END
#        if (sequence_correct[0]==5 and sequence_correct[1]==0) or (sequence_correct[0]==5 and sequence_correct[1]==1):
#            Pe_lower.append(Pe[i,groupStart[1]])
#            St_lower.append(St[i,groupStart[1]])
##        if (sequence_correct[-1]==5 and sequence_correct[-2]==0) or (sequence_correct[-1]==5 and sequence_correct[-2]==1) or (sequence_correct[-1]==0 and sequence_correct[-2]==1):
##            Pe_upper.append(Pe[i,groupEnd[-2]])
##            St_upper.append(St[i,groupEnd[-2]])
        # solution to no solution/unstable
        idx_upper=[]
        a=[1,0]
        b=[2,0]
        c=[3,0]
        d=[4,0]
        e=[1,5]
        f=[2,5]
        g=[3,5]
        h=[4,5]
        for l in range(nSeq):
            if all(sequence_correct[l:l+len(a)] == a)==1:
                idx_upper.append(l+1)
            elif all(sequence_correct[l:l+len(b)] == b)==1:
                idx_upper.append(l+1)
            elif all(sequence_correct[l:l+len(c)] == c)==1:
                idx_upper.append(l+1) 
            elif all(sequence_correct[l:l+len(d)] == d)==1:
                idx_upper.append(l+1)
            elif all(sequence_correct[l:l+len(e)] == e)==1:
                idx_upper.append(l+1) 
            elif all(sequence_correct[l:l+len(f)] == f)==1:
                idx_upper.append(l+1)
            elif all(sequence_correct[l:l+len(g)] == g)==1:
                idx_upper.append(l+1)                   
            elif all(sequence_correct[l:l+len(h)] == h)==1:
                idx_upper.append(l+1)                
        if len(idx_upper)!=0:
#            Pe_upper.append(Pe[i,groupStart[idx_upper[0]]])
#            St_upper.append(St[i,groupStart[idx_upper[0]]])
            Pe_upper.append(Pe[i,groupEnd[idx_upper[0]-1]])
            St_upper.append(St[i,groupEnd[idx_upper[0]-1]]) 
        # no solution to solution
        idx_lower=[]
        a=[5,1]
        b=[5,2]
        c=[5,3]
        d=[5,4]
        e=[5,0]

        for l in range(nSeq):
            if all(sequence_correct[l:l+len(a)] == a)==1:
                idx_lower.append(l+1)
            elif all(sequence_correct[l:l+len(b)] == b)==1:
                idx_lower.append(l+1)
            elif all(sequence_correct[l:l+len(c)] == c)==1:
                idx_lower.append(l+1) 
            elif all(sequence_correct[l:l+len(d)] == d)==1:
                idx_lower.append(l+1)
            elif all(sequence_correct[l:l+len(e)] == e)==1:
                idx_lower.append(l+1) 
        if len(idx_lower)!=0:
#            Pe_upper.append(Pe[i,groupStart[idx_upper[0]]])
#            St_upper.append(St[i,groupStart[idx_upper[0]]])
            Pe_lower.append(Pe[i,groupEnd[idx_lower[0]-1]])
            St_lower.append(St[i,groupEnd[idx_lower[0]-1]]) 
        
        # unstable to stable/partially stable
        idx_unstable=[]    
        a=[0,1] 
        b=[0,2] 
        c=[0,3]
        d=[0,4]
        for l in range(nSeq):
            if all(sequence_correct[l:l+len(a)] == a)==1:
                idx_unstable.append(l+1)
            elif all(sequence_correct[l:l+len(b)] == b)==1:
                idx_unstable.append(l+1)
            elif all(sequence_correct[l:l+len(c)] == c)==1:
                idx_unstable.append(l+1)      
            elif all(sequence_correct[l:l+len(d)] == d)==1:
                idx_unstable.append(l+1)
        if len(idx_unstable)!=0:
#            Pe_unstable.append(Pe[i,groupStart[idx_unstable[0]]])
#            St_unstable.append(St[i,groupStart[idx_unstable[0]]])
            Pe_unstable.append(Pe[i,groupEnd[idx_unstable[0]-1]])
            St_unstable.append(St[i,groupEnd[idx_unstable[0]-1]])
        # Stable to partial
        idx_stable=[]    
        a=[1,3] 
        b=[1,4] 
        c=[2,3]
        d=[2,4]
        for l in range(nSeq):
            if all(sequence_correct[l:l+len(a)] == a)==1:
                idx_stable.append(l+1)
            elif all(sequence_correct[l:l+len(b)] == b)==1:
                idx_stable.append(l+1)
            elif all(sequence_correct[l:l+len(c)] == c)==1:
                idx_stable.append(l+1)                
            elif all(sequence_correct[l:l+len(d)] == d)==1:
                idx_stable.append(l+1)                
        if len(idx_stable)!=0:
            Pe_stable.append(Pe[i,groupStart[idx_stable[0]]])
            St_stable.append(St[i,groupStart[idx_stable[0]]])
        
    # Loop over Stefan
    for j in range(nSt):
        for i in range(nPe-1):
#           No solution to solution 
            if delta[i+1,j]==0 and delta[i,j]==5:
                Pe_unstable_vertical.append(Pe[i,j])
                St_unstable_vertical.append(St[i,j])
#           Unstable to stable - low Le
            if delta[i+1,j]!=0 and delta[i+1,j]!=5 and delta[i,j]==0:
                Pe_stable_vertical.append(Pe[i,j])
                St_stable_vertical.append(St[i,j])                       
    
    return (Pe_stable,St_stable,Pe_unstable,St_unstable,Pe_unstable_vertical,
            St_unstable_vertical,Pe_stable_vertical,St_stable_vertical,
            Pe_upper,St_upper,Pe_lower,St_lower,delta_correct)

    
#%%
    
def getparameterrange(layer_thickness,icb_heatflux,csb_heatflux,thermal_conductivity):
    mol_conc_oxygen_bulk = 8.
    mol_conc_SSi = 8.
    self_diffusion= 0.98e-8
    
    csb_radius=getcsbradius(layer_thickness)
    mass_conc_O,acore=getcsbmassoxygen(mol_conc_oxygen_bulk,mol_conc_SSi)
    deltaV_liquidFeO_solidFe=getchangevolmelting(mol_conc_oxygen_bulk,mol_conc_SSi)
    freezing_speed=getfreezingspeed(icb_heatflux)

    L1,_,_=getL1(csb_radius)
    L2=getL2(mass_conc_O)
    St=getStefan(icb_heatflux,csb_heatflux,csb_radius,freezing_speed)
    Le=getLewis(thermal_conductivity,self_diffusion,deltaV_liquidFeO_solidFe)
    Pe=getPeclet(freezing_speed,csb_radius,self_diffusion,deltaV_liquidFeO_solidFe)

    return L1,L2,St,Le,Pe