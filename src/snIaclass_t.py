#!/usr/bin/env python

import sncosmo
import numpy as np
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catalogs.generation.db import CatalogDBObject
from cStringIO import StringIO
import sys

class SNIa (sncosmo.Model) :

    def __init__(self, galid, snid, ra, dec , z, source = 'salt2'):

        sncosmo.Model.__init__(self, source=source)
        self.snid = snid
        self.galid = galid
        self._ra = ra
        self._dec = dec
        self.z = z

        np.random.seed(self.galid)
    
        c = np.random.normal(0.,0.1)
        hundredyear = 365.0*100.
        x1 = np.random.normal(0.,1.0)
        t0 = np.random.uniform(-hundredyear/2.0, hundredyear/2.0)
        self.set(c=c, x1=x1,  t0=t0, z =z)
        mabs = np.random.normal(-19.3,0.3) 
        self.set_source_peakabsmag(mabs, 'bessellb', 'ab')
        x0 = self.get('x0')

    @property
    def ra(self): 
        return self._ra
    @property
    def dec(self):
        return self._dec


    @property
    def t0(self):
        idx = self.param_names.index("t0")
        return self.parameters[idx]




    



if __name__ == "__main__":


    print sncosmo.__file__
    print SNIa.__class__

    galDB =  CatalogDBObject.from_objid ('galaxyBase')
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
            sn = SNIa( galid=id, snid =id, ra=ra, dec=dec , z=redshift) 

            #print sn.snid, sn.z , sn.ra , sn.dec
            source =  sn._source
            time = 0.0
            phase = (time  - sn.t0) / (1. + sn.z ) 
            print sn.snid, source.flux(phase= time,wave=[ 4000., 5000., 6000.])



