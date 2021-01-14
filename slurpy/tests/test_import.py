#!/usr/bin/env python

"""
Simple test code for slurpy module
"""

from slurpy.coreproperties import icb_radius

def output():
    return icb_radius

def test_output():
    assert output() == 1221500.0
