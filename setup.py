#from ez_setup import use_setuptools
#use_setuptools()
from setuptools import setup, find_packages

setup(name="SNIaCatalogs",
      version="0.1dev",
      description='SNIaCatalogs produces a catalog of simulated SNIa associated to galaxies in a galaxy catalog',
      long_description=''' ''',
      packages=['sniacatalogs','tests', 'examples']
      )
