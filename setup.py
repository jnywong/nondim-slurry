#!/usr/bin/env python

from setuptools import setup

setup(name='slurpy',
      version='1.0',
      description='Python code that solves the 1D, steady, spherical slurry \
      equations outlined in Wong et al (in prep) (see also Wong et al. 2018)',
      url='https://github.com/jnywong/nondim-slurry',
      author='Jenny Wong',
      author_email='jenny.wong@univ-grenoble-alpes.fr',
      license='MIT',
      packages=['slurpy'],
      install_requires=['numpy','scipy','matplotlib','pandas','pickle'],
      # scripts=['slurpy/bin/<filename>.py'],
       package_data={'slurpy': ['lookupdata/*.csv']},
      )
