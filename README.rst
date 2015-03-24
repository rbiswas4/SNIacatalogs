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
  ${HOME}/.local/lib/python2.7/site-packages/mypaths.pth
- LSST Simulation Stack: You will need to use the following branches of the
  current simulation stack.
    - sims_catalogs_generation: feature/ObsMetaData_summary 
    - sims_catalogs_measures: feature/metaData_in_InstanceCatalog
    - sims_photUtils: feature/SIM-1037-only-calculate-requested-magnitudes
- Setup sniacatalogs: From SNIaCatalogs directory:
  >>> python setup.py install --user

TODO: (Major things to be added)
----
These are really part of updates to SNCosmo simulation capabilities

1. Add uncertainties (at some approximation, easy, harder with images)
2. Add intrinsic scatter (will take time)
3. Tune Distributions of SN model parameters

older: (outdated)
-----
src/snIa.py is an instance catalog class to interact with the lsst database
src/snIacLass_t.py is an independent class which interacts with the lsst database
