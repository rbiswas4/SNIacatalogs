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
from snObject import SNObject
import lsst.utils.tests as utilsTests

import utils_for_test as tu
# from unittest import assertAlmostEqual

ra = 204.
dec = -30.
SN = SNObject(ra, dec)
SN.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)

bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
photometry = PhotometryBase()
photometry.loadBandPassesFromFiles(bandPassList)
lsstbands = photometry.bandPassList
l = []
for time in numpy.arange(-20., 50., 1.0):
    t = time*numpy.ones(len(bandPassList))
    t.tolist()
    x = SN.bandMags(lsstbands, time=time)
    # y = SNCosmoSN.bandmag(band=sncosmobands, time=t, magsys='ab')
    e = [time]
    e += x.tolist()
    # e += y.tolist()
    l.append(e)
header = "time(mjd) u g r i z y"
lc = numpy.array(l)
numpy.savetxt('testData/lc.dat', lc, fmt='%10.6f', header=header)

# The 'testData/standard_lc.dat' file was written using the following command
#     for a SN at ra  = 204., -30. , with x0 = 1.846e-6, x1=0.1, c=0., z=0.2
#     using the SALT2-extended model
# numpy.savetxt('testData/standard_lc.dat', numpy.array(l), fmt='%10.6f', header=header)
std_lc = numpy.loadtxt('testData/standard_lc.dat')
numpy.testing.assert_allclose(std_lc, lc)
