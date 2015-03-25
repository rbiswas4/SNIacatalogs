Code to add SNIa to LSST galaxy catalogs:
=========================================

Current Source Files
---------------------
- sniacatalogs/snObject.py : SNObject provides a wrapping of `sncosmo.model`
  adding a few attributes, and methods to use catsim framework to obtain
  fluxes and magnitudes after extinction.
- sniacatalogs/snUniverse.py : Contains rules of distribution of SN model
  parameters.
- sniacatalogs/sncat.py : Instance catalog of SN instantiated using the SNObject class in snObject.py drawing galaxy information from a galaxy catalog.


Usage and Setup:
----------------
- Requires basic SNCosmo (https://github.com/sncosmo/sncosmo)
  sncosmo v1.0.x can be installed via pip. This is also present in the LSST
  simulation stack. The instructions for SNCosmo installation are provided in
  http://sncosmo.readthedocs.org/en/latest/install.html. I typically have a fork
  of the SNcosmo distribution and on the top level SNCosmo directory use:

  >>> python setup.py build 

  and add the path to the created library in build to 
  
  `${HOME}/.local/lib/python2.7/site-packages/mypaths.pth`
- LSST Simulation Stack: You will need to use the following development versions  of packages in the LSST simulation stack. Except in the case of the last two,
  we will be good with the master branch of these package repositories. In the
  last two cases, we would have to use the branches stated here rather than the
  master branch.
    - sims_catalogs_generation  
    - sims_photutils
    - sims_catalogs_measures    
    - sims_utils
    - sims_catutils             
    - sims_coordutils
    - sims_data
    - sims_maf
    - sims_catalogs_generation: feature/ObsMetaData_summary 
    - sims_catalogs_measures: feature/metaData_in_InstanceCatalog
- I am also assuming that you have a anaconda python installed. Some example ipython notebooks can require sqlite3, pandas, seaborn (the last one can be commented out)
- Setup sniacatalogs: Clone the package on your computer, and from SNIaCatalogs directory:

  >>> python setup.py install --user

Example Demonstration
---------------------
After setup, go to the examples directory under the package. Fire up ipython notebook at this location and open the file 
- SNinCatSim.ipynb and evaluate the cells. The notebook should look like http://nbviewer.ipython.org/gist/anonymous/deb1821e8f69def6d7b5 if everything works as expected.
  - This requires sqlite3, pandas, seaborn which are all part of anaconda. If
    Without seaborn, everything will run, so if you do not have seaborn just
    comment out the lines.

  
TODO: (Major things to be added)
================================
These are really part of updates to SNCosmo simulation capabilities

1. Add uncertainties (at some approximation, easy, harder with images)
2. Add intrinsic scatter (will take time)
3. Tune Distributions of SN model parameters

older: (outdated)
-----
src/snIa.py is an instance catalog class to interact with the lsst database
src/snIacLass_t.py is an independent class which interacts with the lsst database
