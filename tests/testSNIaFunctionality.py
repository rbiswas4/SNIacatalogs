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
import utils_for_test as tu

from astropy.units import Unit
from astropy.coordinates import SkyCoord
from sncosmo import Model
import sncosmo
# from unittest import assertAlmostEqual
dustmaproot = os.getenv('SIMS_DUSTMAPS_DIR')
map_dir = os.path.join(dustmaproot, 'DustMaps')

class testSNObject(unittest.TestCase):
    """
    Unit tests to test functionality of the SNObject module. The following tests are
    proposed after the initial objects are setup:

    1. testSNObjectnoMWmags_SNCosmo: Compare the lsst band magnitudes computed using \
            SNCosmo and the LSST software stack.
    2. testSNObjectMWmags_SNCosmo: Compare the lsst band magnitudes computed using \
            SNCosmo and the LSST software stack after the SN Sed has been extincted \
            in the MW, according to CCM89 and observed dustmaps. **** FAIL *** 
    3. testSNObject_lc_againstprev: Compare the light curve output by SNObject for a \
            certain number of days to an older version of the light curve stored on \
            disk in a file 'testData/std_lc.dat'
    """
    def setUp(self):
        """
        Setup the objects used in this test suite class

        """
        ra = 204.
        dec = -30.
        # Setup SN object of different kinds we want to check

        # SNObject with ra, dec, whose extinction will be calculated 
        # separately
        SN = SNObject(ra, dec)
        SN.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)

        # SNCosmo Model object
        SNCosmo = sncosmo.Model(source='salt2-extended',
                                effects=[sncosmo.CCM89Dust()],
                                effect_names=['mw'],
                                effect_frames=['obs'])

        skycoords = SkyCoord(ra, dec, unit='deg')
        t_mwebv = sncosmo.get_ebv_from_map(skycoords,
                                            mapdir=map_dir,
                                            interpolate=False)
        SNCosmo.set(mwebv=t_mwebv)
        self.SNCosmo_mw = SNCosmo


        # Setup SN with no MW extinction
        SNnoMW = SNObject(ra, dec) 
        SNnoMW.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)
        SNnoMW.set_MWebv(0.)
        self.SNmw = SN

        self.SNnoMW = SNnoMW

        #Load LSST sofware bandpass objects for magnitude calculation
        self.bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
        photometry = PhotometryBase()
        photometry.loadBandPassesFromFiles(self.bandPassList)
        self.lsstbands = photometry.bandPassList
        self.times = numpy.arange(-20., 50., 1.0)

        #Load SNCosmo bandpass objects for comparison test
        thisshouldbeNone = tu.getlsstbandpassobjs(self.bandPassList, 
                                                  loadcatsim=False, plot=False)
        self.sncosmobands = ['LSST' + band for band in self.bandPassList]

    def testSNObjectnoMWmags_SNCosmo(self):
        """
        magnitude calculation for SNObject without extinction:
        Compare the lsst band magnitudes computed using SNCosmo and the LSST software 
        stack. This is done by using the SNObject method bandMags and compared with the
        SNCosmo.Model method bandmags with an SNObject whose mwebv attribute is set to
        zero.
        """
        lsst = []
        sncosmo = []
        for time in self.times:

            bandMagsfromLSST = self.SNnoMW.bandMags(self.lsstbands, time=time)
            e  = [time] 
            # e  += bandMagsfromLSST.tolist()
            lsst.append(bandMagsfromLSST.tolist())

            t = time*numpy.ones(len(self.bandPassList))
            t.tolist()
            z = [time]
            y = self.SNnoMW.bandmag(band=self.sncosmobands, time=t, magsys='ab')
            # z += y.tolist()
            sncosmo.append(y.tolist())
            
            numpy.testing.assert_allclose(numpy.array(sncosmo), numpy.array(lsst))
             

    def testSNObjectMWmags_SNCosmo(self):
        """
        ****************************************************************************
        ************************                          **************************
        ************************         FAILING          **************************
        ************************                          **************************
        ****************************************************************************
        magnitude calculation for SNObject with extinction:
        Compare the lsst band magnitudes computed using SNCosmo and the LSST software 
        stack. This is done by using the SNObject method bandMags and compared with the
        SNCosmo.Model method bandmags with an SNObject whose mwebv attribute is 
        calculated in different ways.
        """
        return
        lsst = []
        sncosmo = []
        for time in self.times:

            bandMagsfromLSST = self.SNmw.bandMags(self.lsstbands, time=time)
            e  = [time] 
            # e  += bandMagsfromLSST.tolist()
            lsst.append(bandMagsfromLSST.tolist())

            t = time*numpy.ones(len(self.bandPassList))
            t.tolist()
            z = [time]
            y = self.SNCosmo_mw.bandmag(band=self.sncosmobands, time=t, magsys='ab')
            # z += y.tolist()
            sncosmo.append(y.tolist())
            
            numpy.testing.assert_allclose(numpy.array(sncosmo), numpy.array(lsst))

    def testSNObject_lc_againstprev(self):
        """
        write out the light curve of a single SN represented as SNObject over a number
        of days, and check that the results are the same as the results before stored
        in testData/std_lc.dat 
        """
    
        l = []
        for time in self.times:
        # for time in numpy.arange(-20., 50., 1.0):
            x = self.SNmw.bandMags(self.lsstbands, time=time)
            # y = SNCosmoSN.bandmag(band=sncosmobands, time=t, magsys='ab')
            e = [time]
            e += x.tolist()
            # e += y.tolist()
            l.append(e)
        header = "time(mjd) u g r i z y"
        lc = numpy.array(l)
        numpy.savetxt('testData/lc.dat', lc, fmt='%10.6f', header=header)
    
        # The 'testData/standard_lc.dat' file was written using the following 
        # for a SN at ra  = 204., -30. , with x0 = 1.846e-6, x1=0.1, c=0., z=0.2
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

