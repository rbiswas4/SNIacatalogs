"""
Test the simulation of a SNIa light curve using snObject at specified ra, dec,
and specific SN (SALT2) parameters.
"""
import os
import numpy
import unittest
import matplotlib.pyplot as plt

from lsst.sims.photUtils.Photometry import Sed, Bandpass
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


class testSNObject(unittest.TestCase):
    """
    Unit tests to test functionality of the SNObject module. The following
    tests are proposed after the initial objects are setup:

    1. testSNObjectnoMWmags_SNCosmo: Compare the lsst band magnitudes computed
    using SNCosmo and the LSST software stack when there is no MW extinction.
    2. testSNObjectMWmags_SNCosmo: Compare the lsst band magnitudes computed
    using SNCosmo and the LSST software stack after the SN Sed has been
    extincted in the MW, according to O'Donnell and observed dustmaps.
    3. testSNObject_lc_againstprev: A continuous integration test to compare 
    the light curve output by SNObject for a certain number of days to an 
    older version of the light curve stored on disk in a file 
    'testData/std_lc.dat'
    """
    def setUp(self):
        """
        Setup the objects used in this test suite class. snObject objs are
        required to have a SN sed namely ra, dec, and SALT model parameters
        x0, x1,c and the redshift of the supernova z

        Store objects:

        # NO MW Extinction objects from catsim and SNCosmo
        self.SNnoMW : SN with said parameters, unexincted in MW
        self.SNCosmo_nomw :

        # objects with MW extinction
        self.SNmw : SN with said parameters, exincted in MW
        self.SNCosmo_mw :
        """

        ra = 204.
        dec = -30.
        mjdobs = 571203.
        # Setup SN object of different kinds we want to check

        # In order to compare SN light curves between SNCosmo and Catsim
        # we need to set them up separately

            # Check case with no MW extinction
        SNnoMW = SNObject()
        SNnoMW.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)
        SNnoMW.set_MWebv(0.)
        self.SNnoMW = SNnoMW
        dust = sncosmo.OD94Dust()
        SNCosmo_nomw = sncosmo.Model(source='salt2-extended',
                                     effects=[dust, dust],
                                     effect_names=['host', 'mw'],
                                     effect_frames=['rest', 'obs'])

        SNCosmo_nomw.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)
        SNCosmo_nomw.set(mwebv=0.)
        self.SNCosmo_nomw = SNCosmo_nomw

            ### SNCosmo Model object with MW extinction. Store in self.SNCosmo_mw
        # Object of SNObject class with MW extinction
        SN = SNObject(ra, dec)
        SN.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)
        SNnoMW.set_MWebv(0.)
        self.SNmw = SN

        SNCosmo = sncosmo.Model(source='salt2-extended',
                                effects=[dust, dust],
                                effect_names=['host', 'mw'],
                                effect_frames=['rest', 'obs'])

        SNCosmo.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)
        skycoords = SkyCoord(ra, dec, unit='deg')
        t_mwebv = sncosmo.get_ebv_from_map(skycoords,
                                           mapdir=map_dir,
                                           interpolate=False)
        # Now Set the value of mwebv
        SNCosmo.set(mwebv=t_mwebv)
        self.SNCosmo_mw = SNCosmo

        CCMdust = sncosmo.CCM89Dust()
        SNCosmo_float = sncosmo.Model(source='salt2-extended',
                                      effects=[CCMdust, CCMdust],
                                      effect_names=['host', 'mw'],
                                      effect_frames=['rest', 'obs'])

        SNCosmo_float.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)
        SNCosmo_float.set(mwebv=t_mwebv)
        self.SNCosmo_float = SNCosmo_float
        
            # Load LSST sofware bandpass objects for magnitude calculation
        self.bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
        pbase = PhotometryBase()
        pbase.loadBandpassesFromFiles(self.bandPassList)
        pbase.setupPhiArray_dict()
        self.pbase = pbase
        self.lsstbands = pbase.bandpassDict
        self.times = numpy.arange(-20., 50., 1.0)

        # Load SNCosmo bandpass objects for comparison test
        thisshouldbeNone = tu.getlsstbandpassobjs(self.bandPassList,
                                                  loadcatsim=False,
                                                  plot=False)

        self.sncosmobands = ['LSST' + band for band in self.bandPassList]

    def testSNObjectnoMWmags_SNCosmo(self):
        """
        magnitude calculation for SNObject without extinction:
        Compare the lsst band magnitudes computed using SNCosmo and the LSST
        software stack. This is done by using the SNObject method bandMags and
        compared with the SNCosmo.Model method bandmags with an SNObject whose
        mwebv attribute is set to zero.

        The stored variables are self.SNnomw and self.SNCosmo_nomw
        """
        lsst = []
        sncosmo = []
        for time in self.times:

            bandMagsfromLSST = self.SNnoMW.catsimBandMags(time=time,
                                                                bandpassobject=self.lsstbands,
                                                                phiarray=self.pbase.phiArray).\
                                                                tolist()
            e = [time]
            lsst.append(bandMagsfromLSST)

            t = time*numpy.ones(len(self.bandPassList))
            t.tolist()
            z = [time]

            y = self.SNCosmo_nomw.bandmag(band=self.sncosmobands,
                                          time=t, magsys='ab')
            sncosmo.append(y.tolist())
            numpy.testing.assert_allclose(numpy.array(sncosmo),
                                          numpy.array(lsst))
    def testSNObjectwMWmags_SNCosmo(self):
        """
        magnitude calculation for SNObject without extinction:
        Compare the lsst band magnitudes computed using SNCosmo and the LSST
        software stack. This is done by using the SNObject method bandMags and
        compared with the SNCosmo.Model method bandmags with an SNObject whose
        mwebv attribute is set to zero.

        The stored variables are self.SNnomw and self.SNCosmo_nomw
        """
        lsst = []
        sncosmo = []
        for time in self.times:

            bandMagsfromLSST = self.SNmw.catsimBandMags(time=time,
                                                                bandpassobject=self.lsstbands,
                                                                phiarray=self.pbase.phiArray).\
                                                                tolist()
            e = [time]
            lsst.append(bandMagsfromLSST)

            t = time*numpy.ones(len(self.bandPassList))
            t.tolist()
            z = [time]

            y = self.SNCosmo_mw.bandmag(band=self.sncosmobands,
                                          time=t, magsys='ab')
            sncosmo.append(y.tolist())
            numpy.testing.assert_allclose(numpy.array(sncosmo),
                                          numpy.array(lsst))

    def testSNObject_lc_againstprev(self):
         """
         write out the light curve of a single SN represented as SNObject over a
         number of days, and check that the results are the same as the results
         before stored in testData/std_lc.dat
         """
 
         l = []
         for time in self.times:
             x = self.SNmw.catsimBandMags(time=time, bandpassobject=self.lsstbands,
                                          phiarray=self.pbase.phiArray)

             e = [time]
             e += x.tolist()
             l.append(e)
         header = "time(mjd) u g r i z y"
         lc = numpy.array(l)
         numpy.savetxt('testData/lc.dat', lc, fmt='%10.6f', header=header)
 
         # The 'testData/standard_lc.dat' file was written using the following
         # for a SN at ra  = 204., -30. ,
         # with x0 = 1.846e-6, x1=0.1, c=0., z=0.2
         # using the SALT2-extended model
         # numpy.savetxt('testData/standard_lc.dat', numpy.array(l),
         #               fmt='%10.6f', header=header)
         std_lc = numpy.loadtxt('testData/standard_lc.dat')
         numpy.testing.assert_allclose(std_lc, lc)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(testSNObject)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(shouldExit=True)
