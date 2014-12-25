#!/usr/bin/env python
'''
Example illustrating the writing to DB and reading from database and obtaining light curves of SN using pandas
'''

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
from sncat import SNIaCatalog
from lsst.sims.photUtils.CosmologyObject import CosmologyWrapper 
#import sqliteutils as sq
import testUtilsSNe as sq
import sqlite3
import pandas as pd
wavelenstep = 0.1

import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catalogs.generation.db import ObservationMetaData





def _file2lst(fname, i, mjd):
    d = np.loadtxt(fname, delimiter=',')
    l = list()
    for i, row in enumerate(d):
        obsid = 'obshist' + str(i)
        lst = [obsid] + [mjd] + row.tolist()
        l.append(lst)
    return l
# main example :  create all the 'observations' separated by a day, 


def writeCatalogtoDB(dbfile, dbtable, ascii_root):
    '''
    Write a set of instance catalogs to a sqlite3 database file called db_fname,
    deleting the file if it already exists. This is done in two steps: first the
    write method of instance catalogs is used to write to ASCII files that are
    not removed, and then these are read into the database. The ASCII files 
    are not removed. It is assumed that the directory in which these files 
    are written to exist.
    ''' 
    # erase database if it already exists
    if os.path.exists('testData/sncat.db'):
        print 'deleting previous database'
        os.unlink('testData/sncat.db')
    # Setup connection to write
    connection = sqlite3.connect(dbfile)
    curs = connection.cursor()
    curs.execute('CREATE TABLE if not exists mysncat \
            (id TEXT, mjd FLOAT, snid INT, snra FLOAT, sndec FLOAT,\
            z FLOAT, t0 FLOAT, c FLOAT, x1 FLOAT, x0 FLOAT,\
            mag_u FLOAT, mag_g FLOAT, mag_r FLOAT, mag_i FLOAT,\
            mag_z FLOAT, mag_y FLOAT)')
    # Catalog, and range over which it is written 
    galDB = CatalogDBObject.from_objid('galaxyTiled')
    myMJDS = [570123.15 + 3.*i for i in range(20)]
    for i, myMJD in enumerate(myMJDS):
        myObsMD = ObservationMetaData(boundType='circle',
                                     unrefractedRA=5.0,
                                     unrefractedDec=15.0,
                                     boundLength=0.015,
                                     bandpassName=['u', 'g', 'r', 'i','z', 'y'],
                                     mjd=myMJD)
        catalog = SNIaCatalog(db_obj=galDB,
                              obs_metadata=myObsMD)
        print "====================================="
        print i, type(catalog.usedlsstbands()) , catalog.obs_metadata.mjd
        print "====================================="
        # fname = "data/SNIaCat_" + str(i) + ".txt"
        fname = ascii_root + str(i) + ".txt"
        catalog.write_catalog(fname)
        l = _file2lst(fname, i, mjd=myMJD)
        recs = sq.array2dbrecords(l)
        exec_str = sq.insertfromdata(tablename=dbtable, records=recs,
                                     multiple=True)
        curs.executemany(exec_str, recs)
    connection.commit()
    connection.close()

def getLCsFromDB(dbfile, dbtable, lc_root):
    connection = sqlite3.connect(dbfile)
    df = pd.read_sql('SELECT * FROM ' + dbtable, connection)
    grouped = df.groupby('snid')
    snids = grouped.groups.keys()
    for snid in snids:
        fname =  lc_root + str(snid) + '.dat'
        mylc = df.loc[grouped.groups[snid]]
        mylc.to_csv(fname, na_rep = "NaN", index=False)
    connection.close()
if __name__=="__main__":
    writeCatalogtoDB(dbfile='data/sncat.db', dbtable='mysncat', ascii_root='data/SNIaCat_')
    getLCsFromDB(dbfile='data/sncat.db', dbtable='mysncat', lc_root='data/LightCurves/SN_')




