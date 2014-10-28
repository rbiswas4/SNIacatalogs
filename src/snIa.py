#!\usr/bin/env python

from cStringIO import StringIO
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import DBObject
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

    column_outputs=['snid','z','snra', 'sndec', 'mass_stellar']

    def get_galseed(self) :
        return self.id
    


    def get_snid(self):
        # Not necessarily unique if the same galaxy hosts two SN
        # rethink

        print type(self.column_by_name('id'))
        print self.column_by_name('id')

    
    #@property
    def numobjs(self) :
        return len(self.column_by_name('id'))

    def get_z(self) :
        return self.column_by_name('redshift')

    #dl = cosmo.luminosity_distance(redshift) 
    #da = dl/(1. + redshift**2.0)

    def get_snra(self) :
        numobj = len(self.column_by_name('id'))
        print numobj

        return self.column_by_name('ra') 
        #return np.zeros(numobj)

    def get_sndec(self):
        return self.column_by_name('dec')



    





if __name__=="__main__":

    import lsst.sims.catUtils.baseCatalogModels as bcm
    print bcm.__file__
    from lsst.sims.catalogs.generation.db import ObservationMetaData
    #print ObservationMetaData.__file__
    #circ_bounds = {'ra':5.0, 'dec':15.0,'radius':1.0}
    #print sncosmo.__file__
    #print sncosmo.__version__
    galDB = DBObject.from_objid( 'galaxyBase' )
    #myObsMD = ObservationMetaData(circ_bounds=circ_bounds)
#    myObsMD = ObservationMetaData(boundType='box',
#                                  unrefractedRA=5.0,
#                                  unrefractedDec=15.0,
#                                  boundLength=5.0)
    catalog = SNIaCatalog( db_obj=galDB)
#    catalog = SNIaCatalog( db_obj=galDB, obs_metadata=myObsMD)
    catalog.write_catalog("SNIaCat.txt")

    #catIterator = catalog.query_columns(colnames = ['redshift'])
    #catIterator = galDB.query_columns(colnames = ['redshift'])

    #catalog.write_catalog("testFile.txt")

    




