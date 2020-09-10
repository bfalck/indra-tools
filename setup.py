#!/usr/bin/env python

from distutils.core import setup

setup(name='indratools',
      version='0.4dev6',
      description='Utility code for the Indra Simulations',
      author='Bridget Falck',
      author_email='bridget.falck@jhu.edu',
      packages=['indratools'],
      package_data={'indratools': ['data/*.txt']}
    )
