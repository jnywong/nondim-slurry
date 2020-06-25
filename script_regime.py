#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 17:31:41 2020

@author: jennywong
"""

from slurpy.plot_utils import plot_regime

layer_thickness = 250e3
thermal_conductivity = 100

plot_regime(layer_thickness, thermal_conductivity)