#!/usr/bin/env python

import sncosmo
import numpy as np
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catalogs.generation.db import DBObject
from cStringIO import StringIO
import sys

class SNIa (object) :

    def __init__(self, id, ra, dec, z, mass, modelname):
        self.id = id
        self.ra = ra
        self.dec = dec 
        self._z = z  
        self.mass = mass
        self.modelname = modelname
        
    @property 
    def snid(self) :
        return self.id 

    @property 
    def snra (self) :
        return self.ra 


    @property 
    def sndec (self) :
        return self.dec 

    @property 
    def z (self) :
        return self._z 

    np.random.seed(snid)


    mabs = np.random.normal(-19.3, 0.3)
    if self.modelname == "SALT2exact":
        model = sncosmo.Model(source = "salt2")
    else:
        raise ValueError("Don't know how to treat this supernova modelname")

    model.set(z=z)
    model.set_source_peakabsmag(mabs, 'bessellb', 'ab')
    x0 = model.get('x0')
    p = {'z':z, 't0':uniform(tmin, tmax), 'x0':x0, 'x1': normal(0., 1.), 'c': normal(0., 0.1)}
    params.append(p)


    def dictifySALTparams():
        return 0 

if __name__ == "__main__":


    print sncosmo.__file__
    print SNIa.__class__

    galDB =  DBObject.from_objid ('galaxyBase')
    catalogIterator = galDB.query_columns(colnames=['id','redshift', 'ra', 'dec', 
                                                   'mass_stellar'],
                                          constraint='redshift between 0.01 and 0.1')

    dtype = None
    for chunk in catalogIterator:
        if dtype is None:
            dtype = chunk.dtype
        for record in chunk:
            id = np.int(record['id'])
            mass = np.float(record['mass_stellar'])*10.0**10.0
            ra = np.float(record['ra'])
            dec = np.float(record['dec'])
            redshift = np.float(record['redshift'])
            sn = SNIa( id=id, ra=ra, dec=dec, z=redshift, mass=mass, modelname="SALT2exact") 

            print sn.snid, sn.z , sn.snra , sn.sndec



