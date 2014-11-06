#!\usr/bin/env python
"""
short script to write out an instance catalog from the DB catalog object 
galaxyTiled in a similar way as snIa.py. This is mainly to explore the 
variables and values in the catalog.

R. Biswas
Wed Nov  5 19:33:48 PST 2014
"""

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

    #column_outputs=['snid','snra','sndec', 'z','t0','c', 'x1', 'x0']
    delimiter=' '
    override_formats={'snra':'%8e','sndec':'%8e','c':'%8e','x0':'%8e'}
    column_outputs=['raJ2000','decJ2000',  \
                    'dra', 'ddec',         \
                    'disk_n', 'a_d', 'b_d',\
                    'inc_disk_deg','pa_disk',\
                    'bra', 'bdec',         \
                    'bulge_n', 'a_b', 'b_b',\
                    'inc_bulge_deg', 'pa_bulge',\
                    'u_ab', 'g_ab', 'r_ab', 'i_ab', 'z_ab', 'y_ab']
                    

    def get_snid(self):
        # Not necessarily unique if the same galaxy hosts two SN
        # rethink

        #print type(self.column_by_name('id'))
        return self.column_by_name('id')
    
    @property
    def numobjs(self):
        #return _numobjs
        return len(self.column_by_name('id'))
    
    @property
    def model(self):
        return sncosmo.Model(source='salt2')
    
    #def get_z(self) :
    #    return self.column_by_name('redshift')

    @compound('snra', 'sndec', 'z', 'vra', 'vdec', 'vr')
    def get_angularCoordinates(self) :
        _snra, _sndec, _z =  self.column_by_name('raJ2000'), \
                             self.column_by_name('decJ2000'),\
                             self.column_by_name('redshift')
        #_snra = self.column_by_name('raJ2000') 
        #_sndec = self.column_by_name('decJ2000') 
        #_z = self.column_by_name('redshift')

        #Why does len(_snra) give zero? should it not be of size self.nobj 
        #Yet the addition with np.zeros(self.numobjs) seems to work correctly
        print 'snra', len(_snra)
        _sndec += np.zeros(self.numobjs) 
        _snra +=  np.zeros(self.numobjs) 
        _vra = np.zeros(self.numobjs)
        _vdec = np.zeros(self.numobjs)
        _vr = np.zeros(self.numobjs)
        return ([_snra, _sndec, _z, _vra, _vdec, _vr])

    @compound( 'c', 'x1', 'x0', 't0')
    def get_snparams(self) :

        hundredyear = 100*365.0
        vals  = np.zeros(shape= (self.numobjs,4))
        _z, _id = self.column_by_name('redshift'), self.column_by_name('snid')
        #print _z.min()
        for i,v  in enumerate(vals) :
            np.random.seed(_id[i])
            v[-1 ] = np.random.uniform(-hundredyear/2.0, hundredyear/2.0)
            v[ 0 ] = np.random.normal(0., 0.3 )
            v[ 1 ] = np.random.normal(0., 3.0 ) 
            mabs = np.random.normal( -19.3, 0.3)
            self.model.set(z=_z[i], c=v[0], x1=v[1])
            print 'z val', _z[i]
            self.model.set(z = _z[i])
            print 'from model' , self.model.get('z')
            self.model.set_source_peakabsmag(mabs, 'bessellb', 'ab')
            v[2 ] = v[1]
            v[2 ] = self.model.get('x0')


        return  ([vals[:, 0 ], vals[:, 1], vals[:, 2], vals[:, 3]])
             


    
    
#    def get_c(self):
#        c = np.zeros(self.numobjs, dtype= 'float')
#        for i, id in enumerate(self.column_by_name('id')):
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

#    @property
#    def model(self):
#        return sncosmo.Model(source="SALT2")
#
#    def get_x0(self) :
#        x0 = np.zeros(self.numobjs)
#        for i, id in enumerate(self.column_by_name('id')):
#            np.random.seed(id)
#            mabs = np.random.normal(19.3,0.3) 
#            #z = self.column_by_name('z')[i]
#            #c = self.column_by_name('c')[i]
#            #x1 = self.column_by_name('x1')[i]
#            model = self.model
#            model.set(z=0.5, c=0., x1=0.)
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
    catalog = SNIaCatalog( db_obj=galDB, obs_metadata=myObsMD)
    catalog.write_catalog("galCat.txt")
