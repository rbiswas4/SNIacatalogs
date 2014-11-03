#!\usr/bin/env python

from cStringIO import StringIO
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.generation.db import ObservationMetaData
import sncosmo
from astropy.cosmology import Planck13 as cosmo 

class SNIaCatalog (InstanceCatalog):
    """
    Supernova Type Ia in the catalog are characterized by the  following 
    attributes
    position (ra, dec, redshift), 
    velocity wrt host galaxy, 
    the supernova model (eg. SALT2) 
    and parameters of the supernova model that predict the SED.
    """

    column_outputs=['snid','snra','sndec']
    override_formats={"snra":'%8e',"sndec":'%8e'}

    def get_snid(self):
        # Not necessarily unique if the same galaxy hosts two SN
        # rethink

        #print type(self.column_by_name('id'))
        return self.column_by_name('id')
    
    @property
    def numobjs(self):
        #return _numobjs
        return len(self.column_by_name('id'))
    
    
    def get_z(self) :
        return self.column_by_name('redshift')


    def get_snra(self) :
        _snra = self.column_by_name('raJ2000')
        _snra +=  np.zeros(self.numobjs) 
        return _snra

        

    def get_sndec(self) :
        _sndec = self.column_by_name('decJ2000')
        _sndec += np.zeros(self.numobjs) 
        return _sndec


    
#    def get_c(self):
#        c = np.zeros(self.numobjs, dtype= 'float')
#        for i, id in enumerate(self.column_by_name('id')):
#            np.random.seed(id)
#            c[i] = np.random.normal(0.,0.1)
#        return c
#                
#    def get_x1(self):
#        x1 = np.zeros(self.numobjs, dtype= 'float')
#        for i, id in enumerate(self.column_by_name('id')):
#            np.random.seed(id)
#            x1[i] = np.random.normal(0.,1.0)
#        return x1 
#
#    def get_t0(self) :
#        hundredyear = 100*365.0
#        t0 = np.zeros(self.numobjs, dtype= 'float')
#        for i, id in enumerate(self.column_by_name('id')):
#            np.random.seed(id)
#            t0[i] = np.random.uniform(-hundredyear/2.0, hundredyear/2.0)
#        return t0 
#
#    @property
#    def model(self):
#        return sncosmo.Model(source="SALT2")
#
#    def get_x0(self) :
#        x0 = np.zeros(self.numobjs)
#        for i, id in enumerate(self.column_by_name('id')):
#            np.random.seed(id)
#            mabs = np.random.normal(19.3,0.3) 
#            z = self.column_by_name('z')[i]
#            c = self.column_by_name('c')[i]
#            x1 = self.column_by_name('x1')[i]
#            model = self.model
#            model.set(z=z, c=c, x1=x1)
#            model.set_source_peakabsmag(mabs, 'bessellb', 'ab')
#            x0[i] = -2.5*np.log10(model.get('x0'))
#        return x0
#



if __name__=="__main__":

    import lsst.sims.catUtils.baseCatalogModels as bcm
    print bcm.__file__
    from lsst.sims.catalogs.generation.db import ObservationMetaData
    galDB = CatalogDBObject.from_objid( 'galaxyTiled' )
    myObsMD = ObservationMetaData(boundType='box',
                                  unrefractedRA=5.0,
                                  unrefractedDec=15.0,
                                  boundLength=0.1)
#    catalog = SNIaCatalog( db_obj=galDB)
    catalog = SNIaCatalog( db_obj=galDB, obs_metadata=myObsMD)
    catalog.write_catalog("SNIaCat.txt")

    #catIterator = catalog.query_columns(colnames = ['redshift'])
    #catIterator = galDB.query_columns(colnames = ['redshift'])

    #catalog.write_catalog("testFile.txt")

    




