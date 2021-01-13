#!/usr/bin/env python

import glob
from setuptools import setup

scripts_path = 'scripts/*.py'
# file_list = glob.glob(scripts_path)

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(name='nondim-slurry',
      version='0.0.16',
      description='Python code that solves the 1D, steady, spherical slurry \
      equations outlined in Wong et al (in prep) (see also Wong et al. 2018)',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/jnywong/nondim-slurry',
      license='LICENSE.md',
      author='Jenny Wong',
      author_email='jenny.wong@univ-grenoble-alpes.fr',
      packages=['slurpy'],
      # install_requires=[<pypi_package_name>],
      # scripts=file_list,
      package_data={'slurpy': ['lookupdata/*.csv', 'scripts/*.py']},
      )
