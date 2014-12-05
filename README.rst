Code to add SNIa to LSST galaxy catalogs:
----------------------------------------
Current:
-------
- src/snObject.py : SN Object which has a few more attributes beyond model, and methods to use catsim framework to obtain magnitudes after extinction. 

>>> python snObject.py 
   # should run and produce output light curves for a  SN at z = 0.3 in a file
   # lc.dat.    
   # time, u g r i z y su sg sr si sz sy 
   # u g r i z are calculated using catsim
   # su sg sr si sz sy are calculated using SNCosmo

- src/plotlc.py: plots the light curve output lc.dat comparing them between SNCosmo and the photutils mags. 
- src/sncat.py : Instance catalog of SN instantiated using the SNObject class in snObject.py drawing galaxy information from a galaxy catalog. 

>>> python sncat.py 
# produces instance catalogs for 10 days at an interval of a single day and writes them to disk in files called SNIaCat_i.txt. Each file has the properties of SN observed, and their magnitudes calculated using photUtils. This inlcudes MW extinction through CCM dust relations calculated at the ra, dec positions of the SN.

Usage: 
-----
- Setup LSST using the src/setup.sh provided after modification for your setup. This is different from the usual case,  because we will use a later version of anaconda. The main requirement is due to SNCosmo which uses a later version of astropy for some work. 
- Setup sncosmo: (You only need to do this once). I do this by writing the line: /astro/users/rbiswas/src/sncosmo/build/lib.linux-x86_64-2.7/ to $HOME/.local/lib/python2.7/site-packages/mypaths.pth

TODO: (Major things to be added)
----
1. Add uncertainties (at some approximation, easy)
2. Add intrinsic scatter (will take time) 
3. Add zero point fluctuations 
older: (outdated)
-----
src/snIa.py is an instance catalog class to interact with the lsst database
src/snIacLass_t.py is an independent class which interacts with the lsst database


