#!/usr/bin/env python

from distutils.core import setup
from indratools import __version__

setup(name='indratools',
      version=__version__,
      description='Utility code for the Indra Simulations',
      author='Bridget Falck',
      author_email='bridget.falck@jhu.edu',
      packages=['indratools'],
      package_data={'indratools': ['data/*.txt']}
    )
