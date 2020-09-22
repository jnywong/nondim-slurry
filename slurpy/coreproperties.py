#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 11:35:00 2018

@author: jennywong
"""

###############################################################################
# COREPROPERTIES.PY                                                            #
###############################################################################

# Desciption: estimates of the physical properties of the core

# Time
gigayear = 1e9*365.25*24*60*60; # billion years in seconds
year = 365.25*24*60*60

# PREM
icb_radius = 1.2215e6 # ICB radius
cmb_radius = 3.48e6 # CMB radius
earth_radius = 6.371e6 # radius of the Earth
density_melting = 0.024e4 # change in density upon melting (Alfe 2002, ab initio)
density_solidFe = 1.27e4 # specific density of solid iron
density_liquidFe = density_solidFe-density_melting # specific density of liquid
#                                           iron + S/Si (also reference density)
mass_core=1.94e24 # mass of the core (kg) (Davies, 2015)
mass_oc=1.84e24 # mass of the outer core (kg) (Davies, 2015)
density_liquidO = 0.556e4 # specific density of light element (see Gubbins 2004)
bulk_modulus= 1.324e12 # bulk modulus
hydrostatic_pressure_gradient = 5.509e4 # assumed constant in steady state

# Change in specfic volume (ideal solution theory)
# ...between solid and liquid iron
deltaV_solidFe_liquidFe = 1/density_liquidFe - 1/density_solidFe
# ...between light element and liquid iron
deltaV_liquidO_liquidFe = 1/density_liquidO - 1/density_liquidFe

# Gruneisen parameter
gruneisen=1.51

# Expansion coefficients
alphaT=1e-05 # thermal
alphaXi=1.1 # compositional # CHECK

# Core chemistry
aO=16 # atomic no. of oxygen
aFe=56 # atomic no. of iron
aSSi=30 # atomic no. of silicon/sulphur

# Other properties
heat_capacity=715 # (Jkg^-1K^-1) Gubbins 2003
latent_heat=0.75e6 # (Jkg^-1)
self_diffusion= 0.98e-8
