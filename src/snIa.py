#!\usr/bin/env python

from cStringIO import StringIO
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.generation.db import ObservationMetaData
import sncosmo
from astropy.cosmology import Planck13 as cosmo 

#class SNIaCatalog (object):
class SNIaCatalog (InstanceCatalog):
    """
    Supernova Type Ia in the catalog are characterized by the  following 
    attributes
    position (ra, dec, redshift), 
    velocity wrt host galaxy, 
    the supernova model (eg. SALT2) 
    and parameters of the supernova model that predict the SED.
    """

    column_outputs=['snid','snra','sndec', 'z','t0','c', 'x1', 'x0']
    override_formats={'snra':'%8e','sndec':'%8e','c':'%8e','x0':'%8e'}
    #column_outputs=['raJ2000','decJ2000','snid','z','snra', 'sndec', 'mass_stellar', 'c', 'x1', 't0', "x0"]

    def get_snid(self):
        # Not necessarily unique if the same galaxy hosts two SN
        # rethink

        return self.column_by_name('id')
    
    @property
    def numobjs(self):
        return len(self.column_by_name('id'))
    
    
    @compound('snra', 'sndec', 'z', 'vra', 'vdec', 'vr')
    def get_angularCoordinates(self) :
        _snra, _sndec, _z =  self.column_by_name('raJ2000'), \
                             self.column_by_name('decJ2000'),\
                             self.column_by_name('redshift')
        _sndec += np.zeros(self.numobjs) 
        _snra +=  np.zeros(self.numobjs) 
        _vra = np.zeros(self.numobjs)
        _vdec = np.zeros(self.numobjs)
        _vr = np.zeros(self.numobjs)
        return ([_snra, _sndec, _z, _vra, _vdec, _vr])

    @compound( 'c', 'x1', 'x0', 't0')
    def get_snparams(self) :
        SNmodel = sncosmo.Model(source="salt2")
        hundredyear = 100*365.0
        vals  = np.zeros(shape= (self.numobjs,4))
        _z, _id = self.column_by_name('redshift'), self.column_by_name('snid')
        for i,v  in enumerate(vals) :
            np.random.seed(_id[i])
            v[-1 ] = np.random.uniform(-hundredyear/2.0, hundredyear/2.0)
            v[ 0 ] = np.random.normal(0., 0.3 )
            v[ 1 ] = np.random.normal(0., 3.0 ) 
            mabs = np.random.normal( -19.3, 0.3)
            SNmodel.set(z=_z[i], c=v[0], x1=v[1])
            SNmodel.set_source_peakabsmag(mabs, 'bessellb', 'ab')
            v[2 ] = SNmodel.get('x0')


        return  ([vals[:, 0 ], vals[:, 1], vals[:, 2], vals[:, 3]])

if __name__=="__main__":

    import lsst.sims.catUtils.baseCatalogModels as bcm
    print bcm.__file__
    from lsst.sims.catalogs.generation.db import ObservationMetaData
    galDB = CatalogDBObject.from_objid( 'galaxyTiled' )
    myObsMD = ObservationMetaData(boundType='box',
                                  unrefractedRA=5.0,
                                  unrefractedDec=15.0,
                                  boundLength=0.1)
    catalog = SNIaCatalog( db_obj=galDB, obs_metadata=myObsMD)
    catalog.write_catalog("SNIaCat.txt")

    




