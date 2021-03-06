#!\usr/bin/env python

from cStringIO import StringIO
import sys, os
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
from astropy.cosmology import Planck13 as cosmo



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
    column_outputs = ['snid', 'snra', 'sndec', 'z', 't0', 'c', 'x1', 'x0']
    override_formats = {'snra': '%8e', 'sndec': '%8e', 'c': '%8e', 'x0': '%8e'}
    cannot_be_null = ['x0']
# column_outputs=['raJ2000','decJ2000','snid','z','snra', 'sndec',\
# 'mass_stellar', 'c', 'x1', 't0', "x0"]
    surveyoffset = 570000.0
    SN_thresh  = 100.0

    uband = Bandpass()
    gband = Bandpass()
    rband = Bandpass()
    iband = Bandpass()
    zband = Bandpass()

    lsstbands = [uband, gband, rband, iband, zband]
    banddir = os.path.join(os.getenv('THROUGHPUTS_DIR'), 'baseline')

    bandnames = ['u', 'g', 'r', 'i', 'z']

    def get_mags(self,source, phase ):
        
        sed = Sed()
        print "==========================================="
        print phase
        print "==========================================="
        if phase > -20 and phase < 50 :
            sourceflux = source.flux(phase=phase, wave=self.rband.wavelen*10.)
            sed.setSED(wavelen=self.rband.wavelen, flambda=sourceflux/10.)
        else:
            sed.setSED(wavelen=self.rband.wavelen, flambda=flambda)
        #sed.redshiftSED(redshift=_z[i], dimming=True)
        return [sed.calcMag(bandpass=self.uband),
                sed.calcMag(bandpass=self.gband),
                sed.calcMag(bandpass=self.rband),
                sed.calcMag(bandpass=self.iband),
                sed.calcMag(bandpass=self.zband)]
        # [self.uband, self.gband, self.rband,
        #                 self.iband, self.zband])

    for i, band in enumerate(lsstbands):
        filename = os.path.join(banddir, "total_" + bandnames[i] +".dat")
        band.readThroughput(filename)

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
        SNmodel = sncosmo.Model(source="salt2-extended")
        hundredyear = 100*365.0
        vals = np.zeros(shape=(self.numobjs, 9))
        _z, _id = self.column_by_name('redshift'), self.column_by_name('snid')
        bad = np.nan
        for i, v in enumerate(vals):
            np.random.seed(_id[i])
            v[-1] = np.random.uniform(-hundredyear / 2.0 + self.surveyoffset, hundredyear / 2.0 + self.surveyoffset)
            if np.abs(v[-1] - self.obs_metadata.mjd) > self.SN_thresh:
                #v = np.array([np.nan, np.nan, np.nan, np.nan])
                v[-1] = bad
                v[0] = bad
                v[1] = bad
                v[2] = bad 
                continue
            v[0] = np.random.normal(0., 0.3)
            v[1] = np.random.normal(0., 3.0)
            mabs = np.random.normal(-19.3, 0.3)
            SNmodel.set(z=_z[i], c=v[0], x1=v[1])
            SNmodel.set_source_peakabsmag(mabs, 'bessellb', 'ab')
            v[2] = SNmodel.get('x0')
            phase =  (self.obs_metadata.mjd - v[-1])/(1.0 + _z[i])
            source = SNmodel.source
            v[4:] =  self.get_mags(source, phase = phase)
            # print self.obs_metadata.mjd
        #print self.obs_metadata.bandpass

        return ([vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3],
            vals[:, 4], vals[:,5], vals[:,6], vals[:,7], vals[:,8]])
  
if __name__ == "__main__":

    import lsst.sims.catUtils.baseCatalogModels as bcm
    # import timeit
    print bcm.__file__
    from lsst.sims.catalogs.generation.db import ObservationMetaData
    galDB = CatalogDBObject.from_objid('galaxyTiled')
    myMJDS  = [570123.15 + i for i in range(10)]
    
    for i,myMJD in enumerate(myMJDS):
        myObsMD = ObservationMetaData(boundType='circle',
                                  unrefractedRA=5.0,
                                  unrefractedDec=15.0,
                                  boundLength=0.15,
                                  #boundLength=3.5,
                                  bandpassName=['u','g'],
                                  mjd=myMJD)
        catalog = SNIaCatalog(db_obj=galDB,
                obs_metadata=myObsMD)
        catalog.write_catalog("SNIaCat_" + str(i) + ".txt")
