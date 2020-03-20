#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 10:21:35 2019

@author: wong
"""

###############################################################################
# PLOTREGIME.PY                                                               #
###############################################################################

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy import stats
import os

from scipy.interpolate import interp1d
from matplotlib import cm, colors #, patches

from parametersearch import getphaseboundary, fitphaseboundary
from coreproperties import icb_radius,density_solidFe,latent_heat,deltaV_solidFe_liquidFe
from refparameters import getchangevolmelting
from parameter_range import getdimension

#------------------------------------------------------------------------------
# SELECT INPUTS
#------------------------------------------------------------------------------
 
def getregimedata(L1,L2,Le):
    inputDir='phaseboundaries/sph_L1'+str(np.round(L1,2))+'_L2'+str(np.round(L2,2)) + '_Le'+str(np.round(Le,2))+'/'
    nPe=int(np.genfromtxt(inputDir+'nPe.csv'))
    nSt=int(np.genfromtxt(inputDir+'nSt.csv'))
    Pe=np.genfromtxt(inputDir+'Pe.csv')
    St=np.genfromtxt(inputDir+'St.csv')
    density_jump=np.genfromtxt(inputDir+'density_jump.csv')
    heatflux_cmb=np.genfromtxt(inputDir+'cmb_heatflux.csv')*1e-12
    delta=np.genfromtxt(inputDir+'delta.csv')
    F=np.genfromtxt(inputDir+'F.csv')
    icb_speed=np.genfromtxt(inputDir+'icb_speed.csv')
    
    return (nPe,nSt,Pe,St,density_jump,heatflux_cmb,delta,F,icb_speed)

def smooth_noise(x,y):
    n = len(y)
    x_out = []
    y_out = []
    for i in reversed(range(n)):
        if i==n-1:
            x_out.append(x[i])
            y_out.append(y[i])
            continue
        elif y[i]>y_out[-1]:
            x_out.append(x[i])            
            y_out.append(y[i])
        else:
            while y[i]<y_out[-1]:
                x_out.pop()
                y_out.pop()
            x_out.append(x[i])                
            y_out.append(y[i])            
    
    x_out = x_out[::-1]
    y_out = y_out[::-1]
    
    return x_out,y_out

def fill_gaps(Pe,St,y):
    # Correct numerical artefacts
    nPe = Pe.shape[0]
    ind=[]
    stopOn=0
    for i in range(nPe):
        # Find position of missing data
        for j in range(50,200): # range misses out no solution
            if y[i,j]==0:
                ind.append(j)
        iStart=ind[0]        
        while stopOn==0:
            if y[i,iStart+1]==0:
                ind.append(iStart+1)
                iStart=iStart+1
            else:
                stopOn=1
        # Interpolate data
        xi = np.array([St[i,ind[0]-1],St[i,ind[-1]+1]])
        yi = np.array([y[i,ind[0]-1],y[i,ind[-1]+1]])
        xnew = np.linspace(St[i,ind[0]-1],St[i,ind[-1]+1],len(ind))
        fnew=interp1d(xi,yi)
        for k in range(len(ind)):
            y[i,ind[k]]=fnew(xnew[k])
        ind=[]    
        stopOn=0    
        
        
            
    