#!/usr/bin/env python

from cStringIO import StringIO
import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')

import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.generation.db import ObservationMetaData

from lsst.sims.photUtils import Sed
from lsst.sims.photUtils.Photometry import PhotometryBase as PhotometryBase
from lsst.sims.photUtils.CosmologyObject import CosmologyWrapper 

import astropy
import sncosmo

from snObject import SNObject
from UniversalRules import SNUniverse

import sqlite3

wavelenstep = 0.1

cosmo = CosmologyWrapper()
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
                      'x0'] 
    
    # You can add parameters like fluxes and magnitudes by adding the following
    # variables to the list
    # 'flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y' , 
    # 'mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y']

    override_formats = {'snra': '%8e', 'sndec': '%8e', 'c': '%8e',
            'x0': '%8e'} 
    # You can also change the overrides by 
    # 'flux_u': '%8e', 'flux_g': '%8e', 'flux_r': '%8e',
    # 'flux_i': '%8e', 'flux_z': '%8e', 'flux_y': '%8e'}

    cannot_be_null = ['x0','z', 't0']

    @astropy.utils.lazyproperty
    def mjdobs(self): 
        '''
        The time of observation for the catalog, which is set to be equal
        to obs_metadata.mjd
        '''
        return self.obs_metadata.mjd

    @astropy.utils.lazyproperty
    def badvalues(self):
        '''
        The representation of bad values in this catalog is numpy.nan
        '''
        return np.nan

    @property
    def suppressDimSN(self):
        '''
        Boolean to decide whether to output observations of SN that are too dim
        should be represented in the catalog or not. By default set to True
        '''
        # self._suppressDimSN = True
        if not hasattr(self, '_suppressDimSN'):
            suppressDimSN_default = True
            self._suppressDimSN = suppressDimSN_default
        return self._suppressDimSN

    @suppressDimSN.setter
    def suppressDimSN(self, suppressDimSN):
        '''
        set the value of suppressDimSN of the catalog 

        Parameters
        ----------
        value : Boolean, mandatory
            Value to set suppressDimSN to 
        '''
        # if suppressDimSN is None:
        #    self._suppressDimSN = True
        # else:
        self._suppressDimSN = suppressDimSN
        return self._suppressDimSN
        

    @property
    def suppressHighzSN(self):
        '''
        Boolean to decide whether to output information corresponding to SN 
        at redshift above self.maxz 
        '''

        return True

    @property
    def midSurveyTime(self):
        '''
        The time at the middle of the survey: ie. at the 5 year period.


        .. note: Changing this should not change the statistical
        properties of the survey, but will change the exact SN we find.
        '''
        if not hasattr(self, '_midSurveyTime'):
            midSurveyTime_default = 570000.0
            self._midSurveyTime = midSurveyTime_default
        return self._midSurveyTime


    @midSurveyTime.setter
    def midSurveyTime(self, mymidSurveyTime):
        '''
        set the value of suppressDimSN of the catalog 

        Parameters
        ----------
        value : Boolean, mandatory
            Value to set suppressDimSN to 
        '''
        # if suppressDimSN is None:
        #    self._suppressDimSN = True
        # else:
        self._midSurveyTime = mymidSurveyTime
        return self._midSurveyTime
     
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


        bandPassNames = self.obs_metadata.bandpass
        pbase = PhotometryBase()
        pbase.loadBandpassesFromFiles(bandpassNames=bandPassNames)
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

       
        c, x1, x0, t0, _z , _id, ra, dec = self.column_by_name('c'),\
                                 self.column_by_name('x1'),\
                                 self.column_by_name('x0'),\
                                 self.column_by_name('t0'),\
                                 self.column_by_name('redshift'),\
                                 self.column_by_name('snid'),\
                                 self.column_by_name('raJ2000'),\
                                 self.column_by_name('decJ2000')

        SNobject = SNObject()
        # Initialize return array
        vals = np.zeros(shape=(self.numobjs, 12))
        for i, v in enumerate(vals):
            arr = [_z[i], c[i], x1[i], t0[i], x0[i]]
            testnan = lambda x: x is np.nan
            SNobject.set(z=_z[i], c=c[i], x1=x1[i], t0=t0[i], x0=x0[i]) 
            SNobject.setCoords(ra=ra[i], dec=dec[i])
            SNobject.mwEBVfromMaps()
            # Calculate fluxes
            vals[i, :6] = SNobject.catsimBandFluxes(bandpassobject=self.lsstpbase.bandpassDict,
                                                 time=self.obs_metadata.mjd,
                                            phiarray=self.lsstpbase.phiArray)
            # Calculate magnitudes
            vals[i, 6:] = SNobject.catsimBandMags(bandpassobject=self.lsstpbase.bandpassDict,
                                                time=self.obs_metadata.mjd,
                                            phiarray=self.lsstpbase.phiArray)

        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3],
                vals[:, 4], vals[:, 5], vals[:, 6], vals[:, 7],
                vals[:, 8], vals[:, 9], vals[:, 10], vals[:, 11])
