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
import utils_for_test as tu

from astropy.units import Unit
from astropy.coordinates import SkyCoord
from sncosmo import Model
import sncosmo

import eups
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.generation.db import ObservationMetaData

import testUtilsSNe as sq
import sqlite3
# Global Variables for calculating MWEBV through astropy/SNcosmo
dustmaproot = eups.productDir('sims_dustmaps')
map_dir = os.path.join(dustmaproot, 'DustMaps')


def _file2list(fname, i, mjd):
    d = numpy.loadtxt(fname, delimiter=',')
    l = list()
    for i, row in enumerate(d):
        obsid = 'obshist' + str(i)
        lst = [obsid] + [mjd] + row.tolist()
        l.append(lst)
    return l


class testSNIaCatalog(unittest.TestCase):
    """
    Unit tests to test the functionality of the SNIaCatalog class.

    Requires:
    ---------
        connection to LSST database


    Tests:
    ------
        Write out an instance catalog of SNIa
        Find SNIa 'observed' according to obs_metadata associated\
        with an LSST view, catalog in an instance catalog and write to ascii\
        files testData/SNIaCat_i.txt as output.
        testICatOuput:
        testWriteToSQLite:
        testLCFromSQLite:
    """
    mjds = [570123.15 + 3.*i for i in range(2)]

    @classmethod
    def setUpClass(cls):
        # delete previous test db if present
        if os.path.exists('testData/sncat.db'):
            print 'deleting previous database'
            os.unlink('testData/sncat.db')

        mjds = [570123.15 + 3.*i for i in range(2)]
        galDB = CatalogDBObject.from_objid('galaxyTiled')

        for i, myMJD in enumerate(mjds):
            myObsMD = ObservationMetaData(boundType='circle',
                                          boundLength=0.015,
                                          unrefractedRA=5.0,
                                          unrefractedDec=15.0,
                                          bandpassName=
                                          ['u', 'g', 'r', 'i', 'z', 'y'],
                                          mjd=myMJD)
            catalog = SNIaCatalog(db_obj=galDB, obs_metadata=myObsMD)
            fname = "testData/SNIaCat_" + str(i) + ".txt"
            catalog.write_catalog(fname)

    @classmethod
    def tearDownClass(cls):
        pass

    def testICatOutput(self):
        """
        Check that the output of the instance catalog SNIaCatalog

        """
        stddata = numpy.loadtxt('testData/SNIaCat_0_std.txt', delimiter=',')
        newdata = numpy.loadtxt('testData/SNIaCat_0.txt', delimiter=',')
        stddata.sort(axis=0)
        newdata.sort(axis=0)
        numpy.testing.assert_allclose(stddata, newdata, rtol=1.0e-3)

    def testWriteToSQLite(self):
        """
        """
        connection = sqlite3.connect('testData/sncat.db')
        curs = connection.cursor()
        curs.execute('CREATE TABLE if not exists mysncat (id TEXT, mjd FLOAT,\
                     snid INT, snra FLOAT, sndec FLOAT, z FLOAT, t0 FLOAT,\
                     c FLOAT, x1 FLOAT, x0 FLOAT, flux_u FLOAT, flux_g FLOAT,\
                     flux_r FLOAT, flux_i FLOAT, flux_z FLOAG, flux_y FLOAT,\
                     mag_u FLOAT, mag_g FLOAT, mag_r FLOAT, mag_i FLOAT,\
                     mag_z FLOAT, mag_y FLOAT)')

        for i, myMJD in enumerate(self.mjds):
            fname = "testData/SNIaCat_" + str(i) + ".txt"
            l = _file2list(fname, i, mjd=myMJD)
            recs = sq.array2dbrecords(l)
            exec_str = sq.insertfromdata(tablename='mysncat', records=recs,
                                         multiple=True)
            curs.executemany(exec_str, recs)
        connection.commit()
        curs.execute('SELECT * FROM mysncat')
        print 'LC In Database: '
        lc = curs.fetchall()
        print type(lc)




def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(testSNIaCatalog)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(shouldExit=True)
