#!/usr/bin/env python

from cStringIO import StringIO
import sys
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.generation.db import ObservationMetaData
from lsst.sims.photUtils import Sed
import lsst.sims.photUtils.Bandpass as Bandpass
import sncosmo
from astropy.units import Unit
import astropy.cosmology as cosmology 
#from astropy.cosmology import Planck13 as cosmo
from snObject import SNObject
from lsst.sims.photUtils.CosmologyObject import CosmologyWrapper 
# import sqliteutils as sq
import sqlite3
wavelenstep = 0.1


cosmo = CosmologyWrapper()
# class SNIaCatalog (object):
class SNIaCatalog (InstanceCatalog):
    """
    Supernova Type Ia in the catalog are characterized by the  following
    attributes
    position (ra, dec, redshift),
    velocity wrt host galaxy,
    the supernova model (eg. SALT2)
    and parameters of the supernova model that predict the SED.
    """
    column_outputs = ['snid', 'snra', 'sndec', 'z', 't0', 'c', 'x1',
                      'x0','mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z',
                      'mag_y']
    override_formats = {'snra': '%8e', 'sndec': '%8e', 'c': '%8e',
                      'x0': '%8e'}
    cannot_be_null = ['x0']
# column_outputs=['raJ2000','decJ2000','snid','z','snra', 'sndec',\
# 'mass_stellar', 'c', 'x1', 't0', "x0"]
    surveyoffset = 570000.0
    SN_thresh = 100.0
    def usedlsstbands(self, loadsncosmo=True, loadcatsim=True):

        import eups
        bandPassList = self.obs_metadata.bandpass
        throughputsdir = eups.productDir('throughputs')
        # banddir = os.path.join(os.getenv('THROUGHPUTS_DIR'), 'baseline')
        banddir = os.path.join(throughputsdir, 'baseline')
        lsstbands = []
        lsstbp = {}

        # wavelenstep = 0.1
        for band in bandPassList:
            # setup sncosmo bandpasses
            bandfname = banddir + "/total_" + band + '.dat'
            if loadsncosmo:
                # Usually the next two lines can be merged,
                # but there is an astropy bug currently.
                numpyband = np.loadtxt(bandfname)
                sncosmoband = sncosmo.Bandpass(wave=numpyband[:, 0],
                                               trans=numpyband[:, 1],
                                               wave_unit=Unit('nm'),
                                               name='LSST'+band)
                sncosmo.registry.register(sncosmoband, force=True)
            if loadcatsim:
                # Now load LSST bandpasses for catsim
                lsstbp[band] = Bandpass()
                lsstbp[band].readThroughput(bandfname,
                                            wavelen_step=wavelenstep)
                lsstbands.append(lsstbp[band])
        plot = False
        if plot:
            filterfigs, filterax = plt.subplots()
            for band in bandPassList:
                b = sncosmo.get_bandpass('LSST' + band)
                filterax.plot(b.wave, b.trans, '-k', lw=2.0)
                filterax.set_xlabel(r'$\lambda$, ($\AA$)')
                filterax.set_ylabel(r'transmission')
            plt.show()

        if loadcatsim:
            return lsstbands
        else:
            return None

    def get_snid(self):
        # Not necessarily unique if the same galaxy hosts two SN
        # rethink
        return self.column_by_name('id')

    @property
    def numobjs(self):
        return len(self.column_by_name('id'))

    @compound('snra', 'sndec', 'z', 'vra', 'vdec', 'vr')
    def get_angularCoordinates(self):
        _snra, _sndec, _z = self.column_by_name('raJ2000'), \
            self.column_by_name('decJ2000'), \
            self.column_by_name('redshift')
        _sndec += np.zeros(self.numobjs)
        _snra += np.zeros(self.numobjs)
        _vra = np.zeros(self.numobjs)
        _vdec = np.zeros(self.numobjs)
        _vr = np.zeros(self.numobjs)
        return ([_snra, _sndec, _z, _vra, _vdec, _vr])

    @compound('c', 'x1', 'x0', 't0')
    def get_snparams(self):
        lsstbands = self.usedlsstbands()
        ra, dec = self.column_by_name('raJ2000'),\
            self.column_by_name('decJ2000')
        SNmodel = SNObject()
        hundredyear = 100*365.0
        vals = np.zeros(shape=(self.numobjs, 4))
        _z, _id = self.column_by_name('redshift'), self.column_by_name('snid')
        bad = np.nan
        for i, v in enumerate(vals):
            np.random.seed(_id[i])
            v[-1] = np.random.uniform(-hundredyear / 2.0 +
                                      self.surveyoffset,
                                      hundredyear / 2.0 +
                                      self.surveyoffset)
            if np.abs(v[-1] - self.obs_metadata.mjd) > self.SN_thresh:
                # v = np.array([np.nan, np.nan, np.nan, np.nan])
                v[-1] = bad
                v[0] = bad
                v[1] = bad
                v[2] = bad
                continue
            v[0] = np.random.normal(0., 0.3)
            v[1] = np.random.normal(0., 3.0)
            mabs = np.random.normal(-19.3, 0.3)
            SNmodel.set(z=_z[i], c=v[0], x1=v[1], t0=v[-1])
            # rather than use the SNCosmo function below which uses astropy to calculate
            # distanceModulus, we will use photUtils CosmologyWrapper for consistency
            # SNmodel.set_source_peakabsmag(mabs, 'bessellb', 'ab', cosmo=cosmo)
            mag = mabs + cosmo.cosmology.distanceModulus(_z[i])
            #cosmo.get_cosmologicalDistanceModulus(_z[i])
            SNmodel.source.set_peakmag(mag, band='bessellb', magsys='ab')
            v[2] = SNmodel.get('x0')
            v[3] = v[-1]

        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3]) 

    @compound('mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y')
    def get_snmags(self):
        lsstbands = self.usedlsstbands()
        SNmodel = SNObject()
        vals = np.zeros(shape=(self.numobjs, 6))
        c, x1, x0, t0, _z , _id, ra, dec = self.column_by_name('c'),\
                                 self.column_by_name('x1'),\
                                 self.column_by_name('x0'),\
                                 self.column_by_name('t0'),\
                                 self.column_by_name('redshift'),\
                                 self.column_by_name('snid'),\
                                 self.column_by_name('raJ2000'),\
                                 self.column_by_name('decJ2000')
        for i, v in enumerate(vals):
        #   if x0[i] is np.nan or x0[i] == 'nan':
        #        continue
        #        print ' Got nan'
        #    else:
            SNmodel.set(z=_z[i], c=c[i], x1=x1[i], t0=t0[i], x0=x0[i]) 
            SNmodel.ra=ra[i]
            SNmodel.dec=dec[i]
            SNmodel.mwEBVfromMaps()
            vals[i, :] = SNmodel.bandmags(bandpassobjects=lsstbands,
                                          time=self.obs_metadata.mjd)

        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3], vals[:, 4], vals[:, 5]) 
