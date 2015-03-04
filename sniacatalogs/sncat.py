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

import sqlite3

wavelenstep = 0.1


cosmo = CosmologyWrapper()
class SNIaCatalog (InstanceCatalog, CosmologyWrapper):
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
    surveyoffset = 570000.0
    SN_thresh = 100.0
    maxz = 1.2

    @astropy.utils.lazyproperty
    def mjdobs(self): 
        return self.obs_metadata.mjd

    @astropy.utils.lazyproperty
    def badvalues(self):
        return np.nan

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
    def SN_thresh(self):
        ''' 
        superseded by maxTimeVisibilityforSN
        '''
        return 100.


    @property
    def suppressSNoutsideSurveyTime(self):
        return False

    
    @property
    def maxz(self):
        return 1.2


    @astropy.utils.lazyproperty
    def lsstpbase(self):
        # import eups
        bandPassList = self.obs_metadata.bandpass
        # throughputsdir = eups.productDir('throughputs')
        # banddir = os.path.join(throughputsdir, 'baseline')

        pbase = PhotometryBase()
        pbase.loadBandPassesFromFiles(bandPassList)
        pbase.setupPhiArray_dict()

        return pbase


#    def usedlsstbands(self, loadsncosmo=True, loadcatsim=True):
#
#        import eups
#
#        bandPassList = self.obs_metadata.bandpass
#        throughputsdir = eups.productDir('throughputs')
#        banddir = os.path.join(throughputsdir, 'baseline')
#
#        pbase = PhotometryBase()
#        pbase.loadBandPassesFromFiles(bandPassList)
#        pbase.setupPhiArray_dict()
#
#        lsstbands = []
#        lsstbp = {}
#
#        for band in bandPassList:
#            # setup sncosmo bandpasses
#            bandfname = banddir + "/total_" + band + '.dat'
#            if loadsncosmo:
#                # Usually the next two lines can be merged,
#                # but there is an astropy bug currently.
#                numpyband = np.loadtxt(bandfname)
#                sncosmoband = sncosmo.Bandpass(wave=numpyband[:, 0],
#                                               trans=numpyband[:, 1],
#                                               wave_unit=astropy.units.Unit('nm'),
#                                               name='LSST'+band)
#                sncosmo.registry.register(sncosmoband, force=True)
#
#            if loadcatsim:
#                # Now load LSST bandpasses for catsim
#                lsstbp[band] = Bandpass()
#                lsstbp[band].readThroughput(bandfname,
#                                            wavelen_step=wavelenstep)
#                lsstbands.append(lsstbp[band])
#
#        plot = False
#        if plot:
#            filterfigs, filterax = plt.subplots()
#            for band in bandPassList:
#                b = sncosmo.get_bandpass('LSST' + band)
#                filterax.plot(b.wave, b.trans, '-k', lw=2.0)
#                filterax.set_xlabel(r'$\lambda$, ($\AA$)')
#                filterax.set_ylabel(r'transmission')
#            plt.show()
#
#        if loadcatsim:
#            return lsstbands
#        else:
#            return None

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
        _snra, _sndec, _z = self.column_by_name('raJ2000'), \
            self.column_by_name('decJ2000'), \
            self.column_by_name('redshift')
        _sndec += np.zeros(self.numobjs)
        _snra += np.zeros(self.numobjs)
        _vra = np.zeros(self.numobjs)
        _vdec = np.zeros(self.numobjs)
        _vr = np.zeros(self.numobjs)

        if self.suppressHighzSN:
            _z = np.where(_z > self.maxz, np.nan, _z)
        return ([_snra, _sndec, _z, _vra, _vdec, _vr])

    @compound('c', 'x1', 'x0', 't0')
    def get_snparams(self):
        # RB:  Mon Mar  2 21:02:41 PST 2015
        # Don't see why usedlsstbands is useful for this
        # 
        # lsstbands = self.usedlsstbands()

        SNmodel = SNObject()
        hundredyear = 100*365.0
        vals = np.zeros(shape=(self.numobjs, 4))
        _z, _id, mu  = self.column_by_name('redshift'), self.column_by_name('snid'), self.column_by_name('cosmologicalDistanceModulus')

        bad = np.nan
        for i, v in enumerate(vals):
            np.random.seed(_id[i])
            t0val = np.random.uniform(-hundredyear / 2.0 +
                                      self.midSurveyTime,
                                      hundredyear / 2.0 +
                                      self.midSurveyTime)
            # if np.abs(v[-1] - self.obs_metadata.mjd) > self.SN_thresh:
            #    v[-1] = bad
            #    v[-1] = -1
            #    v[0] = bad
            #    v[1] = bad
            #    v[2] = bad
            #    continue

            # if _z[i] > self.maxz:
            #    v[-1] = bad
            #    v[0] = bad
            #    v[1] = bad
            #    v[2] = bad
            #    continue
            if self.suppressSNoutsideSurveyTime:
                if np.abs(t0val - self.obs_metadata.mjd) > self.SN_thresh:
                    t0val = np.nan

            cval = np.random.normal(0., 0.3)
            x1val = np.random.normal(0., 3.0)
            mabs = np.random.normal(-19.3, 0.3)
            SNmodel.set(z=_z[i], c=cval, x1=x1val, t0=t0val)
            # rather than use the SNCosmo function below which uses astropy to obtain
            # distanceModulus, we will use photUtils CosmologyWrapper for consistency
            # SNmodel.set_source_peakabsmag(mabs, 'bessellb', 'ab', cosmo=cosmo)
            mag = mabs + mu[i]
            SNmodel.source.set_peakmag(mag, band='bessellb', magsys='ab')
            # We can now get x0
            x0val = SNmodel.get('x0')
            vals[i, 0] = cval
            vals[i, 1] = x1val
            vals[i, 2] = x0val
            vals[i, 3] = t0val

        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3]) 

    @compound('flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y','mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y')
    def get_snfluxes(self):

        SNmodel = SNObject()
       
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
            SNmodel.set(z=_z[i], c=c[i], x1=x1[i], t0=t0[i], x0=x0[i]) 
            SNmodel.ra=ra[i]
            SNmodel.dec=dec[i]
            SNmodel.mwEBVfromMaps()
            # Get the `photUtils.SED` object from SNObject
            sed = SNmodel.SNObjectSED(time=self.obs_metadata.mjd, 
                    bandpassobject=self.lsstpbase.bandPassList)
            sed.flambdaTofnu()
            # Calculate fluxes
            vals[i, :6] = sed.manyFluxCalc(phiarray=self.lsstpbase.phiArray,
                    wavelen_step=self.lsstpbase.bandPassList[0].wavelen_step)
            # Calculate magnitudes
            vals[i, 6:] = sed.manyMagCalc(phiarray=self.lsstpbase.phiArray,
                   wavelen_step=self.lsstpbase.bandPassList[0].wavelen_step)

        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3],
                vals[:, 4], vals[:, 5], vals[:, 6], vals[:, 7],
                vals[:, 8], vals[:, 9], vals[:, 10], vals[:, 11])
