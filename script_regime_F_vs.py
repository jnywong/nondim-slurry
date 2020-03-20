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
from scipy import stats
import os

from matplotlib import cm, colors #, patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
from itertools import groupby, count

from parametersearch import getphaseboundary, fitphaseboundary
from parameter_range import getdimension
from regime import getregimedata
from refparameters import getcsbradius, getfreezingspeedfromPe
from coreproperties import self_diffusion

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)
#------------------------------------------------------------------------------
# SELECT INPUTS
#------------------------------------------------------------------------------
# Dimensionless numbers
L1=0.16
L2=0.02
Lewis= [1180.89,354.27]

# Dimensional numbers
layer_thickness= 150e3

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
max_data=8 # contour heat flux - upper limit
min_den = 1e-11
max_den = 1e-10 # contour density jump - upper limit
num_levs=10 # number of contour levels
nTicks=13 # no of ticks on CMB heat flux colorbar
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

fig4,(ax4a,ax4b)=plt.subplots(1,2,figsize=(2*width,height))
fig5,(ax5a,ax5b)=plt.subplots(1,2,figsize=(2*width,height))

for Le in Lewis:
    # Annotation
    fs_ann=10 # fontsize
    if Le==354.27:
        xunstable_ann=50
        xstable_ann=580
        yunstable_ann=2.
        ystable_ann=2
        xpartial_ann=120
        ypartial_ann=4.5
        xnosol_ann=400
        ynosol_ann=0.3
    else:
        xunstable_ann=280
        xstable_ann=800
        yunstable_ann=1.1
        ystable_ann=2
        xpartial_ann=600
        ypartial_ann=4.5
        xnosol_ann=400
        ynosol_ann=0.3
        
    # Polynomial degree for interpolation
    if Le==354.27:
        n_stable1=2
        n_stable2=2
        n_unstable=4
        n_stablelower=1
        n_stableupper=4
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
#    xp_stable,yp_stable=fitphaseboundary(Pe_stable,St_stable,
#                                            Pe_stable[0],Pe_stable[-1],n_stable)
    
    if Le==354.27:
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
    if Le==354.27:
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
               
        
    #%% Contour plot F
    X = Pe
    Y = St
    # Ignore values outside plausible CMB heat flux range
    Z = F
#    Z=ma.masked_where(Y<np.mean(St_stablelower[1:]),Z)
#    Z=ma.masked_where(Y>yUpper+0.8,Z)
#    Z=ma.masked_where(X>xUpper,Z)    
    Z=ma.masked_where(Z<=0,Z)
#    Z = ma.masked_where(heatflux_cmb<5,Z)
#    Z = ma.masked_where(heatflux_cmb>17,Z)    
    print('Max(F) = {:.1f} and min(F) = {:.1f}'.format(Z.compressed().max(),Z.compressed().min()))
    levs = np.linspace(min_data,max_data,num_levs)
    cmap = cm.get_cmap('RdPu',np.size(levs))
    if Le==354.27:
        im3b=ax3b.contourf(X,Y,Z,levs,norm=colors.BoundaryNorm(levs, cmap.N), cmap=cmap)
        x_fill = Pe_unstable
#        ax3b.fill_between(x_fill,0,St_unstable,facecolor='white')
    else:
        im3a=ax3a.contourf(X,Y,Z,levs,norm=colors.BoundaryNorm(levs, cmap.N), cmap=cmap)
#        Pe_unstable.insert(0,m[0])
#        St_unstable.insert(0,np.mean(St_stablelower[1:]))
#        ax3a.fill_between(Pe_unstable,0,St_unstable,facecolor='white') 
        
    #%% StPeL2/L1 vs F
    x = (St*Pe)/Le # + int0/m0 #*L1)/(Le*L2)
    x = ma.masked_where(Y<np.mean(St_stablelower[1:]),x)
    x = ma.masked_where(Pe>xUpper,x)
    x = ma.masked_where(St>yUpper+0.8,x)
    x = ma.masked_where(heatflux_cmb>17,x)
    x = ma.masked_where(heatflux_cmb<5,x)    
    Z = ma.masked_where(heatflux_cmb<5,Z)
    Z = ma.masked_where(heatflux_cmb>17,Z)
    cmap=plt.cm.YlOrRd_r(np.linspace(0.4,1,x.shape[0]))    

    x_fit = np.array([])
    z_fit = np.array([])    

#    i0=60
#    i1=61
    Pe_1d=[]
    St_1d=[]    
    mu=[]
    eta=[]
    k=0
    for i in range(x.shape[0]):   
#    for i in range(i0,i1):
        i2 = np.nonzero(Z[i,:])
        # Ignore cases where there are no solutions        
        if i2[0].size==0:
            continue
        elif x[i,i2].compressed().size!=Z[i,i2].compressed().size:
            continue
#        print(i)
        # Filter out numerical artifacts
        def as_range(g):
            l = list(g)
            return l[0], l[-1]    
        i2 = [as_range(g) for _, g in groupby(i2[0], key=lambda n, c=count(): n-next(c))]
        i2 = np.arange(i2[0][0],i2[0][1],1)              
        ax4a.plot(x[i,i2].compressed(),np.exp(Z[i,i2]).compressed(),color=cmap[i])        
        s_out, i_out = np.polyfit(x[i,i2].compressed(),np.exp(Z[i,i2]).compressed(),1)
        Pe_1d.append(Pe[i,0])
#        St_1d=St[i,:]
        mu.append(s_out)
        eta.append(i_out)
        y_out = s_out*x[i,i2].compressed()+i_out
        ax4a.plot(x[i,i2].compressed(),y_out,'k--')                
        ax4b.plot(x[i,i2].compressed(),Z[i,i2].compressed(),color=cmap[i])
        x_fit =np.append(x_fit,x[i,i2].compressed())
        z_fit =np.append(z_fit,np.exp(Z[i,i2]).compressed())


    idx_sort = np.argsort(x_fit)
    x_fit = x_fit[idx_sort]
    z_fit = z_fit[idx_sort]
    slope, intercept = np.polyfit(x_fit,z_fit,1)
#    print('Slope is {:.2f}, intercept is {:.2f}'.format(slope,intercept))
    y = slope*x_fit+intercept
    ax4a.plot(x_fit,y,'k--')
    y = np.log(y)        
    ax4b.plot(x_fit,y,'k--',label = r'$F=ln({:.2f}x+{:.2f})$'.format(slope,intercept))
    ax4a.set(xlabel=r'$x=\frac{Pe St Li_\xi}{Li_p}$',ylabel=r'$\exp(F)$')
    ax4b.set(xlabel=r'$x=\frac{Pe St Li_\xi}{Li_p}$',ylabel=r'$F$')    
    ax4b.legend()

    print('Lim(F) is {:.2f}'.format(y[-1]))      
        
    #%% Contour plot ICB SPEED
    X = Pe
    Y = St
    # Ignore values outside plausible CMB heat flux range
    csb_radius = getcsbradius(layer_thickness)
    freezing_speed = getfreezingspeedfromPe(Pe,csb_radius,self_diffusion)
    Z = icb_speed/freezing_speed
    # Z=ma.masked_where(Y<np.mean(St_stablelower[1:]),Z)
    Z=ma.masked_where(Y>yUpper+0.8,Z)
    Z=ma.masked_where(X>xUpper,Z)  
    Z=ma.masked_where(Z<=0,Z)
    Z = ma.masked_where(heatflux_cmb<5,Z)
    Z = ma.masked_where(heatflux_cmb>17,Z)
    print('Max(v_s) = {:.3e} and min(v_s) = {:.3e}'.format(Z.compressed().max(),Z.compressed().min()))
    # lev_exp = np.arange(np.floor(np.log10(min_den)-1),
                    # np.ceil(np.log10(max_den))+1)
    # levs = np.power(10,lev_exp)
    cmap = cm.get_cmap('gist_gray',num_levs)
    # Blank out unstable area
    if Le==354.27:
        im3d=ax3d.contourf(X,Y,Z,cmap=cmap) #levs) #, norm=colors.LogNorm()) 
        # ax3d.fill_between(Pe_unstable,0,St_unstable,facecolor='white')
    else:
        im3c=ax3c.contourf(X,Y,Z,cmap=cmap) # levs) #, norm=colors.LogNorm()) 
#        ax3c.fill_between(np.linspace(0,xp_unstable[0]+1,50),np.ones(50)*yLower,np.ones(50)*yUpper,facecolor='white')  
        # ax3c.fill_between(Pe_unstable,0,St_unstable,facecolor='white')

    #%% StPe/Le vs snow speed
    x=St
    Z=icb_speed/freezing_speed
    x = ma.masked_where(Z<=0,x)
    Z = ma.masked_where(Z<=0,Z)
    x = ma.masked_where(Pe>xUpper,x)
    x = ma.masked_where(St>yUpper+0.8,x)
    Z = ma.masked_where(Pe>xUpper,Z)
    Z = ma.masked_where(St>yUpper+0.8,Z)    
    x = ma.masked_where(heatflux_cmb<5,x)
    x = ma.masked_where(heatflux_cmb>17,x)    
    Z = ma.masked_where(heatflux_cmb<5,Z)
    Z = ma.masked_where(heatflux_cmb>17,Z)
    cmap=plt.cm.gist_gray(np.linspace(0.4,1,x.shape[1]))    

    x_fit = np.array([])
    z_fit = np.array([])    

    for i in range(x.shape[1]):    
        i2 = np.nonzero(Z[:,i])
        if i2[0].size==0:
            continue
        i2 = [as_range(g) for _, g in groupby(i2[0], key=lambda n, c=count(): n-next(c))]
        i2 = np.arange(i2[0][0],i2[0][1],1)         
        ax5a.plot(x[i2,i].compressed(),np.exp(Z[i2,i]).compressed(),color=cmap[i])
        ax5b.plot(x[i2,i].compressed(),Z[i2,i].compressed(),color=cmap[i])
        x_fit =np.append(x_fit,x[i2,i].compressed())
        # z_fit =np.append(z_fit,np.exp(Z[i2,i]).compressed())
        z_fit =np.append(z_fit,Z[i2,i].compressed())

    idx_sort = np.argsort(x_fit)
    x_fit = x_fit[idx_sort]
    z_fit = z_fit[idx_sort]
    slope, intercept = np.polyfit(x_fit,z_fit,1)
    y = slope*x_fit+intercept
    # ax5a.plot(x_fit,y,'k--')
    # y = np.log(y)      
    # ax5b.plot(x_fit,y,'k--',label = r'$v_s=ln({:.2e}x+{:.2f})$'.format(slope,intercept))
    ax5b.plot(x_fit,y,'k--',label = r'$v={:.2f}x+{:.2f}$'.format(slope,intercept))    
    ax5a.set(xlabel=r'$x=\frac{St}{Le}$',ylabel=r'$\exp(v)$')
    ax5b.set(xlabel=r'$x=\frac{St}{Le}$',ylabel=r'$v$')    
    ax5a.set_yscale('log')
    # ax5b.set_yscale('log')    
    ax5b.legend()  
    
    print('Lim(v_s) is {:.2e}'.format(y[-1]))      

    if Le==354.27:
        ax3b.text(xnosol_ann,ynosol_ann,r'NO SLURRY',fontsize=fs_ann)
        ax3b.text(xnosol_ann+250,ypartial_ann,r'NO SLURRY',fontsize=fs_ann)
        ax3b.text(xunstable_ann,yunstable_ann,'UNSTABLE',fontsize=fs_ann,rotation=90)
        ax3b.text(xstable_ann,ystable_ann,r'STABLE',fontsize=fs_ann)
        ax3b.text(xpartial_ann,ypartial_ann,'PARTIALLY STABLE',fontsize=fs_ann)
        ax3d.text(xnosol_ann,ynosol_ann,r'NO SLURRY',fontsize=fs_ann)
        ax3d.text(xnosol_ann+250,ypartial_ann,r'NO SLURRY',fontsize=fs_ann)
        ax3d.text(xunstable_ann,yunstable_ann,'UNSTABLE',fontsize=fs_ann,rotation=90)
        ax3d.text(xstable_ann,ystable_ann,r'STABLE',fontsize=fs_ann)
        ax3d.text(xpartial_ann,ypartial_ann,'PARTIALLY STABLE',fontsize=fs_ann)
    else:
        ax3a.text(xnosol_ann,ynosol_ann,r'NO SLURRY',fontsize=fs_ann)
        ax3a.text(xunstable_ann,yunstable_ann,'UNSTABLE',fontsize=fs_ann)
        ax3a.text(xstable_ann,ystable_ann,r'STABLE',fontsize=fs_ann)
        ax3a.text(xpartial_ann,ypartial_ann,'PARTIALLY STABLE',fontsize=fs_ann)
        ax3c.text(xnosol_ann,ynosol_ann,r'NO SLURRY',fontsize=fs_ann)
        ax3c.text(xunstable_ann,yunstable_ann,'UNSTABLE',fontsize=fs_ann)
        ax3c.text(xstable_ann,ystable_ann,r'STABLE',fontsize=fs_ann)
        ax3c.text(xpartial_ann,ypartial_ann,'PARTIALLY STABLE',fontsize=fs_ann)        
        
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
cbar3b.set_label(r'F')
tick_marks= levs #np.linspace(min_den,np.round(max_den,-1),int((max_den-min_den)/20+1))
cbar3d = fig3.colorbar(im3d,cax = cax3d) #, ticks=tick_marks)
cbar3d.set_label(r'ICB speed ($ms^{-1}$)')

fig3.tight_layout(pad=0.3)

#%% SAVE
if saveOn==1:
    saveName='sph_L1'+str(np.round(L1,0)).replace('.0','')+'_L2'+str(np.round(L2,0)).replace('.0','') + '_Le'+str(np.round(Le,0)).replace('.0','')
    if not os.path.exists('phaseboundaries/'+saveName+'/'):
        os.makedirs('phaseboundaries/'+saveName+'/')

    if not os.path.exists('figures/regime_diagram/'):
        os.makedirs('figures/regime_diagram/')
    if rawOn==1:
        fig1.savefig('figures/regime_diagram/'+saveName+'_raw_multi.png',format='png',dpi=200,bbox_inches='tight')
        fig2.savefig('figures/regime_diagram/'+saveName+'_raw_clean_multi.png',format='png',dpi=200,bbox_inches='tight')
    fig3.savefig('figures/regime_diagram/multi_F_vs.png',format='png',dpi=200,bbox_inches='tight') 
    print('Figures saved in figures/regime_diagram/{}'.format(saveName))
    
plt.show()