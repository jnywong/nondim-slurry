#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 12:02:37 2018

@author: jennywong
"""

import numpy as np
from numpy import genfromtxt
from scipy.interpolate import interp1d

folder="slurpy/lookupdata/"

def premgravity(radius):
    # Read PREM data from CSV files
    prem_radius= genfromtxt(folder+'radPREM.csv', delimiter=',')
    prem_gravity=genfromtxt(folder+'gravPREM.csv', delimiter=',')
    
    # Interpolate PREM data
    f=interp1d(prem_radius[14:37],prem_gravity[14:37],kind='cubic',fill_value="extrapolate")
    # start from index 14 to avoid repeating ICB radius at index 13
    
    # Determine CSB density from PREM
    gravity=f(radius)
    
    return gravity

def premdensity(radius):
    # Read PREM data from CSV files
    prem_radius= genfromtxt(folder+'radPREM.csv', delimiter=',')
    prem_density=genfromtxt(folder+'denPREM.csv', delimiter=',')
    
    # Interpolate PREM data
    f=interp1d(prem_radius[14:37],prem_density[14:37],kind='cubic',fill_value="extrapolate")
    # start from index 14 to avoid repeating ICB radius at index 13
    
    # Determine CSB density from PREM
    density=f(radius)
    
    return density
    
def vpspeed(radius):
    # Read PREM data from CSV files
    prem_radius= genfromtxt(folder+'radPREM.csv', delimiter=',')
    prem_vp=genfromtxt(folder+'vpPREM.csv', delimiter=',')
    
    # Interpolate PREM data
    f=interp1d(prem_radius[14:37],prem_vp[14:37],kind='cubic',fill_value="extrapolate")
    # start from index 14 to avoid repeating ICB radius at index 13
    
    # Determine CSB density from PREM
    pwave_speed=f(radius)
    
    return pwave_speed
    
def prempressure(radius): # pressure given radius

    # Read PREM data from CSV files
    prem_radius= genfromtxt(folder+'radPREM.csv', delimiter=',')
    prem_pressure=genfromtxt(folder+'presPREM.csv', delimiter=',')
    
    # Interpolate PREM data
    from scipy.interpolate import interp1d
    f=interp1d(prem_radius[14:37],prem_pressure[14:37],kind='cubic',fill_value="extrapolate")
    # start from index 14 to avoid repeating ICB radius at index 13
    
    # Determine CSB density from PREM
    pressure=f(radius)
    
    return pressure

def premradius(pressure): # radius given pressure
    from numpy import genfromtxt
    # Read PREM data from CSV files
    prem_radius= genfromtxt(folder+'radPREM.csv', delimiter=',')
    prem_pressure=genfromtxt(folder+'presPREM.csv', delimiter=',')
    
    # Interpolate PREM data
    from scipy.interpolate import interp1d
    f=interp1d(prem_pressure[14:37],prem_radius[14:37],kind='cubic',fill_value="extrapolate")
    # start from index 14 to avoid repeating ICB radius at index 13
    
    # Determine CSB density from PREM
    radius=f(pressure)
    
    return radius
    
def liquidus(input_radius): # Obtain liquidus temperature given radius
    # Liquidus polynomial coefficients for pure iron w.r.t. to pressure,
    #                               as in Davies 2015, supplementary material)
    a=np.zeros(4)
    a[0]=1698.55
    a[1]=27.3351
    a[2]=-0.0664736
    a[3]=7.94628e-5

    pressure=prempressure(input_radius)*1e-9 # (GPa)
    
    # Pure iron minus depression from oxygen
    depression=700 # Alfe (2002) EPSL
    tempout=a[0]+a[1]*pressure+a[2]*pressure**2+a[3]*pressure**3-depression
    
    return tempout  

