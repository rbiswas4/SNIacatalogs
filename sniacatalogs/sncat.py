#!/usr/bin/env python

from cStringIO import StringIO
import sys
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.generation.db import ObservationMetaData

from lsst.sims.photUtils import Sed
import lsst.sims.photUtils.Bandpass as Bandpass
from lsst.sims.photUtils.Photometry import PhotometryBase as PhotometryBase
from lsst.sims.photUtils.CosmologyObject import CosmologyWrapper 

import astropy
import sncosmo

from snObject import SNObject
from UniversalRules import SNUniverse

import sqlite3

wavelenstep = 0.1

cosmo = CosmologyWrapper()
class mybandpass(Bandpass, PhotometryBase):
    pass
class SNIaCatalog (InstanceCatalog, CosmologyWrapper, SNUniverse):
    """
    Supernova Type Ia in the catalog are characterized by the  following
    attributes

    Attributes
    ----------
    position : 3-tuple of floats
              (ra, dec, redshift),
    velocity : 3 tuple of floats
              velocity wrt host galaxy in Km/s,
    the supernova model (eg. SALT2)
    and parameters of the supernova model that predict the SED.
    """

    # t_0, c, x_1, x_0 are parameters characterizing a SALT based SN model
    # as defined in sncosmo 
    column_outputs = ['snid', 'snra', 'sndec', 'z', 't0', 'c', 'x1',
                      'x0','flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z',
                      'flux_y' , 'mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z',
                      'mag_y']
    override_formats = {'snra': '%8e', 'sndec': '%8e', 'c': '%8e',
            'x0': '%8e', 'flux_u': '%8e', 'flux_g': '%8e', 'flux_r': '%8e',
            'flux_i': '%8e', 'flux_z': '%8e', 'flux_y': '%8e'}

    cannot_be_null = ['x0','z', 't0']
    # surveyoffset = 570000.0
    # SN_thresh = 100.0
    # maxz = 1.2

    @astropy.utils.lazyproperty
    def mjdobs(self): 
        return self.obs_metadata.mjd

    @astropy.utils.lazyproperty
    def badvalues(self):
        return np.nan

    @property
    def suppressDimSN(self):
        return True

    @property
    def suppressHighzSN(self):
        return True

    @property
    def midSurveyTime(self):
        '''
        The time at the middle of the survey: ie. at the 5 year period.


        .. note: Changing this should not change the statistical
        properties of the survey, but will change the exact SN we find.
        '''
        return 570000.0

    @property
    def maxTimeSNVisible(self):
        '''
        The catalog will provide values for SN flux (even if zero according to
        model) for peak \pm maxTimeVisibilityforSN 
        '''
        return 100.

    
    @property
    def maxz(self):
        return 1.2


    @astropy.utils.lazyproperty
    def lsstpbase(self):
        # import eups
        bandPassNames = self.obs_metadata.bandpass
        # throughputsdir = eups.productDir('throughputs')
        # banddir = os.path.join(throughputsdir, 'baseline')

        # pbase = PhotometryBase(Bandpass)
        pbase = PhotometryBase()
        # pbase = mybandpass()
        pbase.loadBandPassesFromFiles(bandPassNames=bandPassNames)
        pbase.setupPhiArray_dict()

        return pbase


    def get_snid(self):
        # Not necessarily unique if the same galaxy hosts two SN
        return self.column_by_name('id')

    @property
    def numobjs(self):
        return len(self.column_by_name('id'))

    @compound('snra', 'sndec', 'z', 'vra', 'vdec', 'vr')
    def get_angularCoordinates(self):
        '''
        Obtain the coordinates and velocity of the SN from the host galaxy

        Returns
        -------
        `np.ndarray` of coordinara, dec, z, vra, vdec, and vr

        '''
        hostra , hostdec, hostz = self.column_by_name('raJ2000'),\
                                  self.column_by_name('decJ2000'),\
                                  self.column_by_name('redshift')
        snra, sndec, snz, snvra, snvdec, snvr = self.SNCoordinatesFromHost(
                                                    hostra, hostdec, hostz)

        return ([snra, sndec, snz, snvra, sndec, snvr])

    @compound('c', 'x1', 'x0', 't0')
    def get_snparams(self):
        hostz, hostid, hostmu = self.column_by_name('redshift'),\
                                self.column_by_name('snid'),\
                                self.column_by_name('cosmologicalDistanceModulus')

        vals = self.SNparamDistfromHost( hostz, hostid, hostmu)


        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3]) 

    @compound('flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y','mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y')
    def get_snfluxes(self):

        SNobject = SNObject()
       
        c, x1, x0, t0, _z , _id, ra, dec = self.column_by_name('c'),\
                                 self.column_by_name('x1'),\
                                 self.column_by_name('x0'),\
                                 self.column_by_name('t0'),\
                                 self.column_by_name('redshift'),\
                                 self.column_by_name('snid'),\
                                 self.column_by_name('raJ2000'),\
                                 self.column_by_name('decJ2000')

        # Set return array
        vals = np.zeros(shape=(self.numobjs, 12))
        for i, v in enumerate(vals):
            arr = [_z[i], c[i], x1[i], t0[i], x0[i]]
            testnan = lambda x: x is np.nan
            # if any(map(testnan, arr)):
                # vals[i, :] = np.array([np.nan]*5)
                # print 'bad SN'
                # continue
            SNobject.set(z=_z[i], c=c[i], x1=x1[i], t0=t0[i], x0=x0[i]) 
            # SNobject.ra=ra[i]
            # SNobject.dec=dec[i]
            SNobject.setCoords(ra=ra[i], dec=dec[i])
            SNobject.mwEBVfromMaps()
            # Get the `photUtils.SED` object from SNObject
            # sed = SNmodel.SNObjectSED(time=self.obs_metadata.mjd, 
            #        bandpassobject=self.lsstpbase.bandPassList)
            # sed.flambdaTofnu()
            # Calculate fluxes
            vals[i, :6] = SNobject.catsimBandFluxes(bandpassobject=self.lsstpbase.bandPassList,
                                                 time=self.obs_metadata.mjd,
                                            phiarray=self.lsstpbase.phiArray)
            # Calculate magnitudes
            vals[i, 6:] = SNobject.catsimBandMags(bandpassobject=self.lsstpbase.bandPassList,
                                                 time=self.obs_metadata.mjd,
                                            phiarray=self.lsstpbase.phiArray)
            # vals[i, :6] = sed.manyFluxCalc(phiarray=self.lsstpbase.phiArray,
            #        wavelen_step=self.lsstpbase.bandPassList[0].wavelen_step)
            # Calculate magnitudes
            # vals[i, 6:] = sed.manyMagCalc(phiarray=self.lsstpbase.phiArray,
            #      wavelen_step=self.lsstpbase.bandPassList[0].wavelen_step)

        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3],
                vals[:, 4], vals[:, 5], vals[:, 6], vals[:, 7],
                vals[:, 8], vals[:, 9], vals[:, 10], vals[:, 11])
