"""
Test the simulation of a SNIa light curve using snObject at specified ra, dec,
and specific SN (SALT2) parameters. 
"""
import os
import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord

import sncosmo
from sncosmo import Model

from lsst.sims.photUtils.Photometry import PhotometryStars, Sed, Bandpass
from lsst.sims.photUtils.EBV import EBVbase
from lsst.sims.photUtils.Photometry import PhotometryBase
from sniacatalogs.snObject import SNObject

import test_utils as tu
ra = 204.
dec = -30.
SN = SNObject(ra, dec)
SN.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)
print SN
SNCosmoSN = SNObject(ra, dec)
SNCosmoSN.set(x0=1.847e-6, x1=0.1, z=0.2, mwebv=SN._mwebv)

thisshouldbeNone = tu.getlsstbandpassobjs(loadcatsim=False, plot=True)
bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
sncosmobands = ['LSSTu', 'LSSTg', 'LSSTr', 'LSSTi', 'LSSTz', 'LSSTy']
photometry = PhotometryBase()
photometry.loadBandPassesFromFiles(bandPassList)
lsstbands = photometry.bandPassList
print type(lsstbands)
w = sncosmo.get_bandpass(sncosmobands[0]).wave
l = []
for time in np.arange(-20., 50., 1.0):
    t = time*np.ones(len(sncosmobands))
    t.tolist()
    x = SN.lsstbandmags(lsstbands, time=time)
    y = SNCosmoSN.bandmag(band=sncosmobands, time=t, magsys='ab')
    e = [time]
    e += x.tolist()
    e += y.tolist()
    l.append(e)
header = "time(mjd) u g r i z y su sg sr si sz sy"
np.savetxt('test_data/lc.dat', np.array(l), fmt='%10.6f', header=header)
