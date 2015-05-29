"""
Test the simulation of a SNIa light curve using snObject at specified ra, dec,
and specific SN (SALT2) parameters.
"""
import os
import numpy
import unittest
import matplotlib.pyplot as plt

from lsst.sims.photUtils.Photometry import PhotometryStars, Sed, Bandpass
from lsst.sims.photUtils.Photometry import PhotometryBase
import lsst.utils.tests as utilsTests

from snObject import SNObject
from sncat import SNIaCatalog
#import utils_for_test as tu

from astropy.units import Unit
from astropy.coordinates import SkyCoord
from sncosmo import Model
import sncosmo

import eups
from lsst.sims.catalogs.measures.instance import InstanceCatalog
# from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.generation.db import ObservationMetaData

import testUtilsSNe as sq
import sqlite3
# Global Variables for calculating MWEBV through astropy/SNcosmo
dustmaproot = eups.productDir('sims_dustmaps')
map_dir = os.path.join(dustmaproot, 'DustMaps')


mjds = [571203.15]
galDB = CatalogDBObject.from_objid('galaxyTiled')
# galDB = CatalogDBObject.from_objid('galaxyTiled')
class galCopy(InstanceCatalog):
    column_outputs = ['id', 'raJ2000', 'decJ2000', 'redshift']
    override_formats = {'raJ2000': '%8e', 'decJ2000': '%8e'}
for i, myMJD in enumerate(mjds):
    myObsMD = ObservationMetaData(boundType='circle',
                                  boundLength=0.015,
                                  unrefractedRA=5.0,
                                  unrefractedDec=15.0,
                                  bandpassName=['u', 'g', 'r', 'i', 'z', 'y'],
                                  mjd=myMJD)
    gals  = galCopy(db_obj=galDB, obs_metadata=myObsMD)
    gals.write_catalog('gals.dat')
    catalog = SNIaCatalog(db_obj=galDB, obs_metadata=myObsMD)
    print '--- BEFORE ----'
    print '--- AFTER----'
    # print SNIaCatalog.suppressSNoutsideSurveyTime
    print SNIaCatalog.suppressHighzSN
    fname = "SNIaCat_" + str(i) + ".txt"
    print fname, myObsMD.mjd
    catalog.write_catalog(fname)
    SNIaCatalog.suppressHighzSN = True
    fname = "SNIaCat_zsupp" + ".txt"
    print fname, myObsMD.mjd
    catalog.write_catalog(fname)

    SNIaCatalog.suppressSNoutsideSurveyTime = True
    fname = "SNIaCat_zTsupp" + ".txt"
    print fname, myObsMD.mjd
    catalog.write_catalog(fname)
    SNIaCatalog.suppressHighzSN = False
    fname = "SNIaCat_Tsupp" + ".txt"
    print fname, myObsMD.mjd
    catalog.write_catalog(fname)
