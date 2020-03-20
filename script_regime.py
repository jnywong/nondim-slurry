#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:27:42 2020

@author: wong
"""
###############################################################################
# PLOTREGIME.PY                                                               #
###############################################################################

import numpy as np
import numpy.ma as ma
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.font_manager
from scipy import stats
import os
import pickle

from matplotlib import cm, colors #, patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
from itertools import groupby, count
from scipy.interpolate import interp1d

from refparameters import getcsbradius, getfreezingspeedfromPe
from parametersearch import getphaseboundary, fitphaseboundary
from parameter_range import getdimension, getdimensionless
from regime import getregimedata, smooth_noise
from coreproperties import self_diffusion

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)
#------------------------------------------------------------------------------
# SELECT INPUTS
#------------------------------------------------------------------------------
# Dimensionless numbers
thermal_conductivity= [100,30] #,30]
layer_thickness= 150e3
Lewis=[]
for n in thermal_conductivity:
    (L1,L2,_,_,lewis) = getdimensionless(layer_thickness,1,1,n)
    Lewis.append(lewis)

# Plotting parameters
aspectRatio=0.75 # figure aspect ratio
rawOn=0 # plot raw data
contour_hf=0 # contour heat flux
contour_density=0 # contour density jump
dimensionalOn=0 # add dimensional contours

# limits
xLower=0
xUpper=1000
yLower=0
yUpper=5

min_data=5 # contour heat flux - lower limit
max_data=17 # contour heat flux - upper limit
min_den = 60
max_den = 340 # contour density jump - upper limit
min_v = 0 #1e-12 # ICB speed
max_v = 20 # 1e-10
num_levs=100 # number of contour levels
nTicks=13 # no of ticks on CMB heat flux colorbar
den_ticks = 40 # step size for density jump colorbar
n_stableupperLower=14

# Max Pe range for max density jump
Pemax = xUpper

# Save figures?
saveOn=1

#------------------------------------------------------------------------------

width, height = plt.figaspect(aspectRatio)
#fig3,((ax3a,ax3b),(ax3c,ax3d))=plt.subplots(2,2,figsize=(2*width,2.5*height),sharex=True,sharey=True)
fig3 = plt.figure(figsize=(1.75*width,2*height))
ax3a = fig3.add_subplot(2,2,1)
ax3b = fig3.add_subplot(2,2,2,sharex=ax3a,sharey=ax3a)
ax3c = fig3.add_subplot(2,2,3,sharex=ax3a,sharey=ax3a)
ax3d = fig3.add_subplot(2,2,4,sharex=ax3a,sharey=ax3a)

fig4,((ax4a,ax4b),(ax4c,ax4d))=plt.subplots(2,2,figsize=(2*width,2*height)) #,sharex = 'col')
fig5,(ax5a,ax5b)=plt.subplots(1,2,figsize=(2*width,height))

# Format scientific notation to 10^x
f_sci = mticker.ScalarFormatter(useOffset=False, useMathText=True)
g_sci = lambda x,pos : "${}$".format(f_sci._formatSciNotation('%1.2e' % x))
fmt = mticker.FuncFormatter(g_sci)

for Le in Lewis:
    # Annotation
    fs_ann=10 # fontsize
    if Le==Lewis[1]:
        xunstable_ann=50
        xstable_ann=300
        yunstable_ann=2.
        ystable_ann=2
        xpartial_ann=500
        ypartial_ann=4.5
        xnosol_ann=400
        ynosol_ann=0.3
    else:
        xunstable_ann=280
        xstable_ann=600
        yunstable_ann=1.1
        ystable_ann=3
        xpartial_ann=480
        ypartial_ann=4.5
        xnosol_ann=400
        ynosol_ann=0.3
        
    # Polynomial degree for interpolation
    if Le==Lewis[1]:
        n_stable1=2
        n_stable2=2
        n_unstable=2
        n_stablelower=2
        n_stableupper=2
    else:
        n_stable=4
        n_unstable=2
        n_stablelower=1
        n_stableupper=2        
    
    # %% Read data
    (nPe,nSt,Pe,St,density_jump,heatflux_cmb,delta,F,icb_speed)= \
                                                    getregimedata(L1,L2,Le)
    
    # %% Determine phase boundaries of the regime diagram
    (Pe_stable,St_stable,Pe_unstable,St_unstable,Pe_unstable_vert,St_unstable_vert,
         Pe_stable_vert,St_stable_vert,Pe_stableupper,St_stableupper,Pe_stablelower,
               St_stablelower,delta_correct) = getphaseboundary(Pe,St,delta,nPe,nSt)
     
    # %%Plot raw data
    if rawOn==1:
        fig1,ax1=plt.subplots(1,1,figsize=(width,height))
        ax1.scatter(Pe,St,c=delta)
        ax1.set_xlim((xLower,xUpper))
        ax1.set_ylim((yLower,yUpper))
        
        fig2=plt.figure
        fig2,ax2=plt.subplots(1,1,figsize=(width,height))
        ax2.scatter(Pe,St,c=delta_correct)
        ax2.set_xlim((xLower,xUpper))
        ax2.set_ylim((yLower,yUpper))
        
        ax3a.plot(Pe_stablelower,St_stablelower,'rx',linewidth=2.0,solid_capstyle='round')
        ax3a.plot(Pe_stableupper,St_stableupper,'rx',linewidth=2.0,solid_capstyle='round')
        ax3a.plot(Pe_stable,St_stable,'rx',linewidth=2.0,solid_capstyle='round')
        ax3a.plot(Pe_unstable,St_unstable,'rx',linewidth=2.0,solid_capstyle='round')
        ax3a.plot(Pe_unstable_vert,St_unstable_vert,'rx',linewidth=2.0,solid_capstyle='round')
    
    # %% Interpolate phase boundaries
   
    if Le==Lewis[1]:
        xp_stableupper,yp_stableupper=fitphaseboundary(Pe_stableupper[45:],St_stableupper[45:],
                  Pe_stableupper[0],xUpper,n_stableupper)
        xp_stablelower,yp_stablelower=fitphaseboundary(Pe_stablelower[55:65],St_stablelower[55:65],
                  Pe_stablelower[55]+5,Pe_stablelower[-1],n_stablelower)          
        xp_stable1,yp_stable1=fitphaseboundary(Pe_stable[:7],St_stable[:7],
                                        Pe_stable[0]-5,Pe_stable[7],n_stable1)            
        xp_stable2,yp_stable2=fitphaseboundary(Pe_stable[7:58],St_stable[7:58],
                                        Pe_stable[7],xUpper,n_stable2)  
        # Append lower to unstable and stablelower
        St_lower = St_stablelower[0]            
        idx = int(np.argwhere(Pe_stablelower>Pe_unstable[-1])[0])
        Pe_unstable.append(Pe_stablelower[idx])
        St_unstable.append(St_lower)    
        np.insert(xp_stablelower,0,Pe_stablelower[55])
        np.insert(yp_stablelower,0,St_lower)
        # Clean data
#        Pe_stableupper,St_stableupper = smooth_noise(Pe_stableupper[33:64],St_stableupper[33:64])
#        del St_stableupper[3:7]
#        del Pe_stableupper[3:7]     
#        del St_stableupper[0]
#        del Pe_stableupper[0]              
        xp_unstable = Pe_unstable[4:10]
        yp_unstable = St_unstable[4:10]
        xp_unstable.append(xp_unstable[-1]+2)
        yp_unstable.append(St_lower)
        xp_unstable.insert(0,0)
        yp_unstable.insert(0,yUpper)
#        ax3b.plot(Pe_stable,St_stable,'rx',linewidth=2.0,solid_capstyle='round')
    else:
        # Append first point of lower to upper
        Pe_stableupper.insert(0,Pe_unstable[0])
        St_stableupper.insert(0,St_unstable[0])    
        xp_stableupper,yp_stableupper=fitphaseboundary(Pe_stableupper[51:72],St_stableupper[51:72],
                    Pe_stableupper[50],Pe_stableupper[-1],n_stableupper)
        xp_unstable,yp_unstable=fitphaseboundary(Pe_unstable[13:], \
                    St_unstable[13:],Pe_unstable[13],Pe_unstable[-1], n_unstable)
        xp_stable1,yp_stable1=fitphaseboundary(Pe_stable,St_stable,
                                            Pe_stable[0]-10,Pe_stable[-1],n_stable)
    
    # %% Plot phase boundaries 
    
    if Le==Lewis[1]:
        ax3b.plot(xp_stableupper,yp_stableupper,'k',linewidth=2.0,solid_capstyle='round',
                 linestyle='-')
        ax3b.plot(xp_stablelower,yp_stablelower,'k',linewidth=2.0,solid_capstyle='round',
                 linestyle='-')        
        ax3b.plot(xp_unstable,yp_unstable,'k',linewidth=2.0,solid_capstyle='round')
        ax3b.plot(xp_stable1,yp_stable1,'k',linewidth=2.0,solid_capstyle='round',linestyle='--')
        ax3b.plot(xp_stable2,yp_stable2,'k',linewidth=2.0,solid_capstyle='round',linestyle='--')
        ax3b.hlines(St_lower,xmin=0,xmax = Pe_stablelower[55]+5,linewidth=2.0,
                    color='k',linestyle='-')  
        ax3d.plot(xp_stableupper,yp_stableupper,'k',linewidth=2.0,solid_capstyle='round',
                 linestyle='-')
        ax3d.plot(xp_stablelower,yp_stablelower,'k',linewidth=2.0,solid_capstyle='round',
                 linestyle='-')           
        ax3d.plot(xp_unstable,yp_unstable,'k',linewidth=2.0,solid_capstyle='round',linestyle='-')
        ax3d.plot(xp_stable1,yp_stable1,'k',linewidth=2.0,solid_capstyle='round',linestyle='--')
        ax3d.plot(xp_stable2,yp_stable2,'k',linewidth=2.0,solid_capstyle='round',linestyle='--')        
        ax3d.hlines(St_lower,xmin=0,xmax = Pe_stablelower[55]+5,linewidth=2.0,
                    color='k',linestyle='-')            
    else:
        m,_=stats.mode(Pe_unstable_vert[1:])
        ax3a.plot(Pe_unstable,St_unstable,'k',linewidth=2.0,solid_capstyle='round')
        ax3a.vlines(m,np.mean(St_stablelower[1:]),yUpper,linewidth=2.0,color='k',linestyle='-')    
        ax3a.plot(xp_stable1,yp_stable1,'k',linewidth=2.0,solid_capstyle='round',linestyle='--')
        ax3a.plot(xp_stableupper,yp_stableupper,'k',linewidth=2.0,solid_capstyle='round',
                 linestyle='-')
        ax3c.plot(xp_stable1,yp_stable1,'k',linewidth=2.0,solid_capstyle='round',linestyle='--')
        ax3c.plot(Pe_unstable,St_unstable,'k',linewidth=2.0,solid_capstyle='round')
        ax3c.plot(xp_stableupper,yp_stableupper,'k',linewidth=2.0,solid_capstyle='round',
                 linestyle='-')
        ax3c.vlines(m,np.mean(St_stablelower[1:]),yUpper,linewidth=2.0,color='k',linestyle='-')  
        ax3c.hlines(np.mean(St_stablelower[1:]),m,xUpper,
                   linewidth=2.0,color='k',linestyle='-')
        ax3a.hlines(np.mean(St_stablelower[1:]),m,xUpper,
                   linewidth=2.0,color='k',linestyle='-')        
            
    #%% Contour plot CMB HEAT FLUX
#    if contour_hf==1:    
    X = Pe
    Y = St
    # Ignore values outside plausible CMB heat flux range
    Z = np.nan_to_num(heatflux_cmb, neginf = max_data)
    Z[Z<0]=max_data
    Z[Z>max_data] = max_data
#    Z=ma.masked_less(Z,min_data)
#    Z=ma.masked_greater(Z,max_data)
#    Z=ma.masked_where(delta==5,Z)
#    Z=ma.masked_where(delta==0,Z)        
    levs = np.linspace(min_data,max_data,num_levs)
    cmap = cm.get_cmap('YlOrRd_r',np.size(levs))
    if Le==Lewis[1]:
        Pe_unstable.insert(0,0)
        St_unstable.insert(0,0)        
        x_fill = Pe_unstable
        y_fill = St_unstable
#        ax3b.fill_between(xp_stableupper,np.mean(St_stablelower),yp_stableupper,facecolor='floralwhite')
#        ax3b.fill_between(x_fill,y_fill,yUpper,facecolor='floralwhite')
#        ax3b.fill_between(np.linspace(x_fill[-1],xp_stableupper[0]),St_lower,yUpper,facecolor='floralwhite')        
        im3b=ax3b.contourf(X,Y,Z,levs,norm=colors.BoundaryNorm(levs, cmap.N),cmap=cmap,extend='both')        
        ax3b.fill_between(xp_unstable,St_lower,yp_unstable,facecolor='gainsboro')  
        x_fill = np.linspace(xLower,xUpper)
        ax3b.fill_between(x_fill,0,St_lower,facecolor='white')
        x_fill = xp_stableupper
        ax3b.fill_between(x_fill,yp_stableupper,yUpper,facecolor='white')
        x_fill = xp_stablelower
        ax3b.fill_between(x_fill,yLower,yp_stablelower,facecolor='white')
   
    else:
        Pe_unstable.insert(0,m[0])
        St_unstable.insert(0,np.mean(St_stablelower[1:]))
#        ax3a.fill_between(xp_stableupper,np.mean(St_stablelower),yp_stableupper,facecolor='floralwhite')        
#        ax3a.fill_between(Pe_unstable,St_unstable,yUpper,facecolor='floralwhite')
        x_fill = np.linspace(xp_unstable[-1],xp_stableupper[0])
#        ax3a.fill_between(x_fill,np.mean(St_stablelower),yUpper,facecolor='floralwhite')        
        im3a=ax3a.contourf(X,Y,Z,levs,norm=colors.BoundaryNorm(levs, cmap.N),cmap=cmap,extend='both')  
        ax3a.fill_between(Pe_unstable,np.mean(St_stablelower),St_unstable,facecolor='gainsboro') 
        ax3a.fill_between(np.linspace(m[0],Pe_unstable[1]),np.mean(St_stablelower),yUpper,facecolor='gainsboro')        
        x_fill = np.linspace(0,xUpper)
        ax3a.fill_between(x_fill,0,np.mean(St_stablelower),facecolor='white')
        x_fill = np.linspace(0,m[0])        
        ax3a.fill_between(x_fill,0,yUpper,facecolor='white')        
        ax3a.fill_between(xp_stableupper,yp_stableupper,yUpper,facecolor='white') 
        
        
    #%% x vs CMB heat flux
    x = (St*Pe)/(Le)
    x = ma.masked_where(Y<np.mean(St_stablelower[1:]),x)
    x = ma.masked_where(Pe>xUpper,x)
    x = ma.masked_where(St>yUpper+0.8,x)
    Z=ma.masked_where(Y<np.mean(St_stablelower[1:]),Z)
    Z=ma.masked_where(Y>yUpper+0.8,Z)
    Z=ma.masked_where(X>xUpper,Z)    
    Z=ma.masked_where(Z<=0,Z)
    Z=ma.masked_where(delta==0,Z)
    x = ma.masked_where(heatflux_cmb>17,x)
    x = ma.masked_where(heatflux_cmb<5,x)    
    Z = ma.masked_where(heatflux_cmb<5,Z)
    Z = ma.masked_where(heatflux_cmb>17,Z)

    x_fit = np.array([])
    z_fit = np.array([])    

#    i0=60
#    i1=61
    for i in range(x.shape[0]):   
#    for i in range(i0,i1):
        i2 = np.nonzero(Z[i,:])
        # Ignore cases where there are no solutions        
        if i2[0].size==0:
            continue
        elif x[i,i2].compressed().size!=Z[i,i2].compressed().size:
            continue
        # Filter out numerical artifacts
        def as_range(g):
            l = list(g)
            return l[0], l[-1]    
        i2 = [as_range(g) for _, g in groupby(i2[0], key=lambda n, c=count(): n-next(c))]
        i2 = np.arange(i2[0][0],i2[0][1],1)        
        if Le == Lewis[0]:
            ax4a.plot(x[i,i2].compressed(),Z[i,i2].compressed(),color='orange')
            ax5a.plot(x[i,i2].compressed(),Z[i,i2].compressed(),color='orange')            
        else:
            ax4a.plot(x[i,i2].compressed(),Z[i,i2].compressed(),color='gold')
            ax5a.plot(x[i,i2].compressed(),Z[i,i2].compressed(),color='gold')            
        x_fit =np.append(x_fit,x[i,i2].compressed())
        z_fit =np.append(z_fit,Z[i,i2].compressed())

    idx_sort = np.argsort(x_fit)
    x_fit = x_fit[idx_sort]
    z_fit = z_fit[idx_sort]
    slope, intercept = np.polyfit(x_fit,z_fit,1)
    y = slope*x_fit+intercept
    if Le == Lewis[0]:    
        ax4a.plot(x_fit,y,'k',ls='--',lw=3,label = r'$Q^c=${}$x{:.2f}$'.format(fmt(slope),intercept))
        ax5a.plot(x_fit,y,'k',ls='--',lw=3,label = r'$Q^c=${}$x{:.2f}$'.format(fmt(slope),intercept))        
    else:
        ax4a.plot(x_fit,y,'k',ls=':',lw=3,label = r'$Q^c=${}$x{:.2f}$'.format(fmt(slope),intercept))
        ax5a.plot(x_fit,y,'k',ls=':',lw=3,label = r'$Q^c=${}$x{:.2f}$'.format(fmt(slope),intercept)) 
        
    # Save scaling
    if saveOn==1:
        if not os.path.exists('scalings/'):
            os.makedirs('scalings/')    
        if Le==Lewis[0]:
            with open('scalings/highLe_Qcmb.pkl', 'wb') as f:
                pickle.dump([slope, intercept],f)
        else:
            with open('scalings/lowLe_Qcmb.pkl', 'wb') as f:
                pickle.dump([slope, intercept],f)           
    
    #%% Contour plot DENSITY JUMP
    X = Pe
    Y = St
    
    # Ignore values outside plausible CMB heat flux range
    Z = density_jump
#    Z=ma.masked_less_equal(Z,0)    
#    Z=ma.masked_greater(Z,max_den)
#    Z=ma.masked_where(heatflux_cmb>max_data,Z)
#    Z=ma.masked_where(heatflux_cmb<min_data,Z)
#    n_levs = int((max_den-min_den)/20+1)
    n_levs = num_levs
    levs = np.linspace(min_den,max_den,n_levs)
#    levs = np.arange(min_den,max_den+den_ticks,den_ticks)
    cmap = cm.get_cmap('YlGnBu_r',np.size(levs))
    # Blank out unstable area
    if Le==Lewis[1]:
        x_fill = Pe_unstable
        y_fill = St_unstable
#        ax3d.fill_between(xp_stableupper,np.mean(St_stablelower),yp_stableupper,facecolor='azure')
#        ax3d.fill_between(x_fill,y_fill,yUpper,facecolor='azure')
#        ax3d.fill_between(np.linspace(x_fill[-1],xp_stableupper[0]),np.mean(St_stablelower),yUpper,facecolor='azure')             
        im3d=ax3d.contourf(X,Y,Z,levs,norm=colors.BoundaryNorm(levs, cmap.N), cmap=cmap) #get,extend='both')
        ax3d.fill_between(xp_unstable,St_lower,yp_unstable,facecolor='gainsboro')  
        x_fill = np.linspace(xLower,xUpper)
        ax3d.fill_between(x_fill,0,St_lower,facecolor='white')
        x_fill = xp_stableupper
        ax3d.fill_between(x_fill,yp_stableupper,yUpper,facecolor='white')        
        x_fill = xp_stablelower
        ax3d.fill_between(x_fill,yLower,yp_stablelower,facecolor='white')        
    else:
#        ax3c.fill_between(xp_stableupper,np.mean(St_stablelower),yp_stableupper,facecolor='azure')        
#        ax3c.fill_between(Pe_unstable,St_unstable,yUpper,facecolor='azure')
        x_fill = np.linspace(xp_unstable[-1],xp_stableupper[0])
#        ax3c.fill_between(x_fill,np.mean(St_stablelower),yUpper,facecolor='azure')              
        im3c=ax3c.contourf(X,Y,Z,levs,norm=colors.BoundaryNorm(levs, cmap.N), cmap=cmap, extend='both')
        ax3c.fill_between(Pe_unstable,np.mean(St_stablelower),St_unstable,facecolor='gainsboro') 
        ax3c.fill_between(np.linspace(m[0],Pe_unstable[1]),np.mean(St_stablelower),yUpper,facecolor='gainsboro')        
        x_fill = np.linspace(0,xUpper)
        ax3c.fill_between(x_fill,0,np.mean(St_stablelower),facecolor='white')
        x_fill = np.linspace(0,m[0])        
        ax3c.fill_between(x_fill,0,yUpper,facecolor='white')        
        ax3c.fill_between(xp_stableupper,yp_stableupper,yUpper,facecolor='white') 

    #%% x vs density jump
    x = Pe/Le
    x = ma.masked_where(Y<np.mean(St_stablelower[1:]),x)
    x = ma.masked_where(Pe>xUpper,x)
    x = ma.masked_where(St>yUpper+0.8,x)
    Z=ma.masked_where(Y<np.mean(St_stablelower[1:]),Z)
    Z=ma.masked_where(delta==0,Z)
    Z=ma.masked_where(Y>yUpper+0.8,Z)
    Z=ma.masked_where(X>xUpper,Z)    
    Z=ma.masked_where(Z<=0,Z)    
#    Z=ma.masked_where(delta==5,Z)
#    Z=ma.masked_where(delta==0,Z)        
    x = ma.masked_where(heatflux_cmb>17,x)
    x = ma.masked_where(heatflux_cmb<5,x)    
    Z = ma.masked_where(heatflux_cmb<5,Z)
    Z = ma.masked_where(heatflux_cmb>17,Z)   

    x_fit = np.array([])
    z_fit = np.array([])    

#    i0=55
#    i1=56
    for i in range(x.shape[1]):   
#    for i in range(i0,i1):
        i2 = np.nonzero(Z[:,i])
        if i2[0].size==0:
            continue
        elif x[i2,i].compressed().size!=Z[i2,i].compressed().size:
            continue
        # Filter out numerical artifacts
        i2 = [as_range(g) for _, g in groupby(i2[0], key=lambda n, c=count(): n-next(c))]
        i2 = np.arange(i2[0][0],i2[0][1],1)
        if Le == Lewis[0]:  
            ax4c.plot(x[i2,i].compressed(),Z[i2,i].compressed(),color='steelblue')
            ax5b.plot(x[i2,i].compressed(),Z[i2,i].compressed(),color='steelblue')            
        else:
            ax4c.plot(x[i2,i].compressed(),Z[i2,i].compressed(),color='turquoise')
            ax5b.plot(x[i2,i].compressed(),Z[i2,i].compressed(),color='turquoise')            
        x_fit =np.append(x_fit,x[i2,i].compressed())
        z_fit =np.append(z_fit,Z[i2,i].compressed())

    idx_sort = np.argsort(x_fit)
    x_fit = x_fit[idx_sort]
    z_fit = z_fit[idx_sort]
    slope, intercept = np.polyfit(x_fit,z_fit,1)
    y = slope*x_fit+intercept
    if Le == Lewis[0]:  
        ax4c.plot(x_fit,y,'k',ls='--',lw=3,label = r'$\Delta \rho={:.2f}x+{:.2f}$'.format(slope,intercept))
        ax5b.plot(x_fit,y,'k',ls='--',lw=3,label = r'$\Delta \rho={:.2f}x+{:.2f}$'.format(slope,intercept))        
    else:
        ax4c.plot(x_fit,y,'k',ls=':',lw=3,label = r'$\Delta \rho={:.2f}x+{:.2f}$'.format(slope,intercept))
        ax5b.plot(x_fit,y,'k',ls=':',lw=3,label = r'$\Delta \rho={:.2f}x+{:.2f}$'.format(slope,intercept))         
    
    # Save scaling
    if saveOn==1:
        if not os.path.exists('scalings/'):
            os.makedirs('scalings/')    
        if Le==Lewis[0]:
            with open('scalings/highLe_den.pkl', 'wb') as f:
                pickle.dump([slope, intercept],f)
        else:
            with open('scalings/lowLe_den.pkl', 'wb') as f:
                pickle.dump([slope, intercept],f)         
    
    #%% Find max density jump 
    # Pe index range
    x0 = np.argwhere(Pe<Pe_unstable[-1])[0][0]
    x1 = np.argwhere(Pe < xUpper)[-1][0] 
    n_max = x1-x0+1
    max_densityjump = np.zeros(n_max)
    max_Pe = np.zeros(n_max)
    max_St = np.zeros(n_max)
    test = np.zeros(n_max)
    k=0
    for i in range(x0,x1+1): 
        # find St within stable regime
#        valid_idx = np.argwhere(St[i,:]< f_stableupper(Pe[i,:]))
        valid_idx = np.argwhere(St[i,:]<yUpper)      
        y0=0
        y1=valid_idx[-1][0]
            
        max_densityjump[k]=Z[i,y0:y1].max()
        if  np.isnan(max_densityjump[k]):
            max_densityjump[k]=0
        test[k]=Z[i,y0:y1].argmax()
        max_Pe[k]=Pe[i,Z[i,y0:y1].argmax()]
        max_St[k]=St[i,Z[i,y0:y1].argmax()]
        k=k+1
    i_max = np.argmax(max_densityjump)    
    if Le==Lewis[1]:
        thermal_conductivity=30
    else:
        thermal_conductivity=100
    max_Qi, max_Qsl = getdimension(layer_thickness,thermal_conductivity,
                                   max_Pe[i_max],max_St[i_max])
    print('Maximum density jump is {:.2f} kg/m^3 at (Pe,St)=({:2f},{:2f}) or (Qi,Qsl)=({:.2f},{:.2f})'.format(
            max_densityjump.max(),max_Pe[i_max],max_St[i_max],max_Qi,max_Qsl))
#    ax3a.plot(max_Pe[i_max],max_St[i_max], '*', color='yellow', markeredgecolor='k', markersize=15)
#    ax3a.plot(Pe_star,St_star, '*', color='yellow', markeredgecolor='k', markersize=15)    
    
    #%% x vs F
    Z=F
    x = ((St*Pe)/Le) #*L2)/(L1)
    x = ma.masked_where(Y<np.mean(St_stablelower[1:]),x)
    x = ma.masked_where(Pe>xUpper,x)
    x = ma.masked_where(St>yUpper+0.8,x)
    x = ma.masked_where(heatflux_cmb>17,x)
    x = ma.masked_where(heatflux_cmb<5,x)    
    Z = ma.masked_where(heatflux_cmb<5,Z)
    Z = ma.masked_where(heatflux_cmb>17,Z)
#    Z=ma.masked_where(delta==5,Z)
#    Z=ma.masked_where(delta==0,Z)        
    cmap=plt.cm.YlOrRd_r(np.linspace(0.4,1,x.shape[0]))    

    x_fit = np.array([])
    z_fit = np.array([])    

#    i0=60
#    i1=61
    for i in range(x.shape[0]):   
#    for i in range(i0,i1):
        i2 = np.nonzero(Z[i,:])
        # Ignore cases where there are no solutions        
        if i2[0].size==0:
            continue
        elif x[i,i2].compressed().size!=Z[i,i2].compressed().size:
            continue
        # Filter out numerical artifacts
        def as_range(g):
            l = list(g)
            return l[0], l[-1]    
        i2 = [as_range(g) for _, g in groupby(i2[0], key=lambda n, c=count(): n-next(c))]
        i2 = np.arange(i2[0][0],i2[0][1],1)        
        if Le == Lewis[0]:
            ax4b.plot(x[i,i2].compressed(),Z[i,i2].compressed(),color='purple')
        else:
            ax4b.plot(x[i,i2].compressed(),Z[i,i2].compressed(),color='orchid')
        x_fit =np.append(x_fit,x[i,i2].compressed())
        z_fit =np.append(z_fit,np.exp(Z[i,i2]).compressed())

    idx_sort = np.argsort(x_fit)
    x_fit = x_fit[idx_sort]
    z_fit = z_fit[idx_sort]
    fit_out,residual,_,_,_ = np.polyfit(x_fit,z_fit,1,full=True)
    slope = fit_out[0]
    intercept = fit_out[1]
#    print('F residual is {:.2f}'.format(residual[0]))
#    print('Slope is {:.2f}, intercept is {:.2f}'.format(slope,intercept))
    y = slope*x_fit+intercept
#    ax4c.plot(x_fit,y,'k--')
    y = np.log(y)        
    if Le == Lewis[0]:    
        ax4b.plot(x_fit,y,'k',ls='--',lw=3,label = r'$F=ln({:.2f}x{:.2f})$'.format(slope,intercept))
    else:
        ax4b.plot(x_fit,y,'k',ls=':',lw=3,label = r'$F=ln({:.2f}x{:.2f})$'.format(slope,intercept))
#    ax4b.set_ylim([y.min(),y.max()])
    ax4b.set_yticks(np.arange(4.5,7.5,0.5))    
    ax4b.set(ylabel=r'$F$',xlabel=r'$x=\frac{Pe St}{Le}$')
    ax4b.legend(loc='lower right', fontsize = 14)  

    print('Lim(F) is {:.2f}'.format(y[-1]))     
    
    # Save scaling
    if saveOn==1:
        if not os.path.exists('scalings/'):
            os.makedirs('scalings/')    
        if Le==Lewis[0]:
            with open('scalings/highLe_F.pkl', 'wb') as f:
                pickle.dump([slope, intercept],f)
        else:
            with open('scalings/lowLe_F.pkl', 'wb') as f:
                pickle.dump([slope, intercept],f)         
    
    #%% x vs ICB speed
    csb_radius = getcsbradius(layer_thickness)
    freezing_speed = getfreezingspeedfromPe(Pe,csb_radius,self_diffusion)
    Z = icb_speed/freezing_speed
    x = St
#    x = ma.masked_where(St<np.mean(St_stablelower[1:]),x)
#    x = ma.masked_where(Pe>xUpper,x) 
#    x = ma.masked_where(St>yUpper+0.8,x)
    x = ma.masked_where(heatflux_cmb>17,x)
    x = ma.masked_where(heatflux_cmb<5,x)  
    Z = ma.masked_where(heatflux_cmb<5,Z)
    Z = ma.masked_where(heatflux_cmb>17,Z)
#    Z=ma.masked_where(delta==5,Z)
#    Z=ma.masked_where(delta==0,Z)        

    x_fit = np.array([])
    z_fit = np.array([])    

    for i in range(x.shape[0]):    
        i2 = np.nonzero(Z[i,:])
        if i2[0].size==0:
            continue
        i2 = [as_range(g) for _, g in groupby(i2[0], key=lambda n, c=count(): n-next(c))]
        i2 = np.arange(i2[0][0],i2[0][1],1)         
#        ax5a.plot(x[i,i2].compressed(),np.exp(Z[i,i2]).compressed(),color=cmap[i])
        if Le == Lewis[0]:
            ax4d.plot(x[i,i2].compressed(),Z[i,i2].compressed(),color='green')
        else:
            ax4d.plot(x[i,i2].compressed(),Z[i,i2].compressed(),color='limegreen')
        x_fit =np.append(x_fit,x[i,i2].compressed())
        # z_fit =np.append(z_fit,np.exp(Z[i,i2]).compressed())
        z_fit =np.append(z_fit,Z[i,i2].compressed())        

    idx_sort = np.argsort(x_fit)
    x_fit = x_fit[idx_sort]
    z_fit = z_fit[idx_sort]
    slope, intercept = np.polyfit(x_fit,z_fit,1)
    y = slope*x_fit+intercept
#    ax5a.plot(x_fit,y,'k--')
    # y = np.log(y)      
    if Le == Lewis[0]:
        ax4d.plot(x_fit,y,'k',ls='--',lw=3,label = r'$v=${}$x$+{:.2f}'.format(fmt(slope),intercept))
    else:
        ax4d.plot(x_fit,y,'k',ls=':',lw=3,label = r'$v=${}$x$+{:.2f}'.format(fmt(slope),intercept))
    ax4d.set(xlabel=r'$x=St$',ylabel=r'$v \ (\mathrm{ms}^{-1})$')    
    # ax4d.set_yscale('log')
    ax4d.set_xlim(0, 5)    
    ax4d.set_ylim(min_v, max_v)
    ax4d.legend(loc='lower right', fontsize = 14)  
    
    # Save scaling
    if saveOn==1:
        if not os.path.exists('scalings/'):
            os.makedirs('scalings/')    
        if Le==Lewis[0]:
            with open('scalings/highLe_v.pkl', 'wb') as f:
                pickle.dump([slope, intercept],f)
        else:
            with open('scalings/lowLe_v.pkl', 'wb') as f:
                pickle.dump([slope, intercept],f)          

#%% Annotate
    if Le==Lewis[1]:
        ax3b.text(xnosol_ann,ynosol_ann,r'NO SLURRY',fontsize=fs_ann)
        ax3b.text(xnosol_ann+400,ypartial_ann,'NO \n SLURRY',fontsize=fs_ann)
        ax3b.text(xunstable_ann,yunstable_ann,'UNSTABLE',fontsize=fs_ann,rotation=90)
        ax3b.text(xstable_ann,ystable_ann,r'STABLE',fontsize=fs_ann)
        ax3b.text(xpartial_ann,ypartial_ann,'PARTIALLY \n STABLE',fontsize=fs_ann)
        ax3d.text(xnosol_ann,ynosol_ann,r'NO SLURRY',fontsize=fs_ann)
        ax3d.text(xnosol_ann+400,ypartial_ann,'NO \n SLURRY',fontsize=fs_ann)
        ax3d.text(xunstable_ann,yunstable_ann,'UNSTABLE',fontsize=fs_ann,rotation=90)
        ax3d.text(xstable_ann,ystable_ann,r'STABLE',fontsize=fs_ann)
        ax3d.text(xpartial_ann,ypartial_ann,'PARTIALLY \n STABLE',fontsize=fs_ann)
    else:
        ax3a.text(xnosol_ann,ynosol_ann,r'NO SLURRY',fontsize=fs_ann)
        ax3a.text(xunstable_ann,yunstable_ann,'UNSTABLE',fontsize=fs_ann)
        ax3a.text(xstable_ann,ystable_ann,r'STABLE',fontsize=fs_ann)
        ax3a.text(xpartial_ann,ypartial_ann,'PARTIALLY \n STABLE',fontsize=fs_ann)
        ax3c.text(xnosol_ann,ynosol_ann,r'NO SLURRY',fontsize=fs_ann)
        ax3c.text(xunstable_ann,yunstable_ann,'UNSTABLE',fontsize=fs_ann)
        ax3c.text(xstable_ann,ystable_ann,r'STABLE',fontsize=fs_ann)
        ax3c.text(xpartial_ann,ypartial_ann,'PARTIALLY \n STABLE',fontsize=fs_ann)        
        
#%% Plot finishing touches
ax3a.set_xlim((xLower,xUpper))
ax3a.set_ylim((yLower,yUpper))

ax3c.set(xlabel=r'$Pe$',ylabel=r'$St$')
ax3a.set(ylabel=r'$St$')
ax3d.set(xlabel=r'$Pe$')

plt.setp(ax3a.get_xticklabels(), visible=False)
plt.setp(ax3b.get_xticklabels(), visible=False)
plt.setp(ax3b.get_yticklabels(), visible=False)
plt.setp(ax3d.get_yticklabels(), visible=False)

# Colorbars
divider3a = make_axes_locatable(ax3a)
cax3a = divider3a.append_axes("right", size="5%", pad=0.05)

divider3b = make_axes_locatable(ax3b)
cax3b = divider3b.append_axes("right", size="5%", pad=0.05)

divider3c = make_axes_locatable(ax3c)
cax3c = divider3c.append_axes("right", "5%", pad=0.05)

divider3d = make_axes_locatable(ax3d)
cax3d = divider3d.append_axes("right", "5%", pad=0.05)

cbar3a = fig3.colorbar(im3a, cax = cax3a)
fig3.delaxes(ax=cax3a)

cbar3c = fig3.colorbar(im3c, cax = cax3c)
fig3.delaxes(ax=cax3c)

tick_marks = np.arange(min_data,max_data+1,2)
cbar3b = fig3.colorbar(im3b, ticks=tick_marks, cax = cax3b)
cbar3b.set_label(r'CMB heat flux ($TW$)')
tick_marks=np.arange(min_den,max_den+den_ticks,den_ticks)
cbar3d = fig3.colorbar(im3d, ticks=tick_marks,cax = cax3d)
cbar3d.set_label(r'density jump ($kg m^{-3}$)')

ax4a.set(ylabel=r'$Q^c \ (\mathrm{TW})$',xlabel=r'$x=\frac{Pe St}{Le}$')
ax4a.set_ylim([min_data,max_data])
ax4a.set_yticks(np.arange(min_data,max_data+1,2))
ax4a.legend(loc='lower right', fontsize = 14)
ax4b.set_xticks(np.arange(0,12,2))
ax4d.set_xticks(np.arange(0,6,1))
ax5a.set_xticks(np.arange(0,12,2))
ax5a.set(ylabel=r'$Q^c \ (\mathrm{TW})$',xlabel=r'$x=\frac{Pe St}{Le}$')
ax5a.set_xticks(np.arange(0,12,2))
ax5a.set_ylim([min_data,max_data])
ax5a.set_yticks(np.arange(min_data,max_data+1,2))
ax5a.legend(loc='lower right', fontsize = 14)    

ax4c.set(xlabel=r'$x=\frac{Pe}{Le}$',ylabel=r'$\Delta \rho \ (\mathrm{kgm}^{-3})$')    
ax4c.set_xlim([0,8])  
ax4c.set_ylim([min_den,max_den])    
ax4c.set_yticks(np.arange(min_den,max_den+20,den_ticks))
ax4c.legend(loc='lower right', fontsize = 14)    
ax5b.set(xlabel=r'$x=\frac{Pe}{Le}$',ylabel=r'$\Delta \rho \ (\mathrm{kgm}^{-3})$')    
ax5b.set_xlim([0,8])
ax5b.set_xticks(np.arange(0,9,2))      
ax5b.set_ylim([min_den,max_den])    
ax5b.set_yticks(np.arange(min_den,max_den+20,den_ticks))
ax5b.legend(loc='lower right', fontsize = 14)       

fig3.tight_layout(pad=0.3)
fig4.tight_layout(pad=0.3)
fig5.tight_layout(pad=0.3)

# Subfigure titles
ax3a.set_title('(a)',x=0.05,y=1,fontsize=14)
ax3b.set_title('(b)',x=0.05,y=1,fontsize=14)
ax3c.set_title('(c)',x=0.05,y=1,fontsize=14)
ax3d.set_title('(d)',x=0.05,y=1,fontsize=14)
ax4a.set_title('(a)',x=0.05,y=1,fontsize=14)
ax4b.set_title('(b)',x=0.05,y=1,fontsize=14)
ax4c.set_title('(c)',x=0.05,y=1,fontsize=14)
ax4d.set_title('(d)',x=0.05,y=1,fontsize=14)

# Anti aliasing for vector images
for c in im3a.collections:
    c.set_edgecolor("face")
for c in im3b.collections:
    c.set_edgecolor("face")
for c in im3c.collections:
    c.set_edgecolor("face")    
for c in im3d.collections:
    c.set_edgecolor("face")
#%% SAVE
if saveOn==1:
    saveName='sph_L1'+str(np.round(L1,0)).replace('.0','')+'_L2'+str(np.round(L2,0)).replace('.0','') + '_Le'+str(np.round(Le,0)).replace('.0','')
#    if not os.path.exists('phaseboundaries/'+saveName+'/'):
#        os.makedirs('phaseboundaries/'+saveName+'/')

    if not os.path.exists('figures/regime_diagram/'):
        os.makedirs('figures/regime_diagram/')
    if rawOn==1:
        fig1.savefig('figures/regime_diagram/'+saveName+'_raw_multi.pdf',format='pdf',dpi=200,bbox_inches='tight')
        fig2.savefig('figures/regime_diagram/'+saveName+'_raw_clean_multi.pdf',format='pdf',dpi=200,bbox_inches='tight')
    fig3.savefig('figures/regime_diagram/multi_Qc_den_v2.pdf',format='pdf',dpi=200,bbox_inches='tight')
    fig4.savefig('figures/regime_diagram/multi_scalings.pdf',format='pdf',dpi=200,bbox_inches='tight')
    fig5.savefig('figures/regime_diagram/Qc_den_scalings.pdf',format='pdf',dpi=200,bbox_inches='tight')    
    print('Figures saved in figures/regime_diagram/{}'.format(saveName))
    
plt.show()