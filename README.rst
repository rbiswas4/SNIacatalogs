Code to add SNIa to LSST galaxy catalogs:
----------------------------------------
Current:
src/sncat.py : Instance catalog
src/snObject.py : SN Object which has a few more attributes, and methods

older:
src/snIa.py is an instance catalog class to interact with the lsst database
src/snIacLass_t.py is an independent class which interacts with the lsst database

Usage: 
-----
- Setup LSST using the src/setup.sh provided after modification for your setup. This is different from the usual case,  because we will use a later version of anaconda. The main requirement is due to SNCosmo which uses a later version of astropy for some work. 
- Setup sncosmo: (You only need to do this once). I do this by writing the line: /astro/users/rbiswas/src/sncosmo/build/lib.linux-x86_64-2.7/ to $HOME/.local/lib/python2.7/site-packages/mypaths.pth

 
