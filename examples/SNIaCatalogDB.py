#!/isr/bin/env python
'''
Example illustrating the following tasks:
    1. Write instance catalogs corresponding to pointings to ASCII files
    2. Read the ASCII files and ingest to a sqlite3 DB
    3. Read the database and group the observations according to SN name to
       create light curves which are now written to disk.
    The example assumes that the directory ${PWD}/data and
    ${PWD}/data/Lightcurves exist
'''

from cStringIO import StringIO
import sys
import os
import numpy as np
import sqlite3
import pandas as pd
# import sncosmo
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
from astropy.units import Unit
import astropy.cosmology as cosmology
from snObject import SNObject
from sncat import SNIaCatalog
from lsst.sims.photUtils.CosmologyObject import CosmologyWrapper
import lsst.sims.catUtils.utils.testUtilsSNe as sq


wavelenstep = 0.1


def _file2lst(fname, i, mjd):
    """
    creates a lst of observations from a Catalog written out
    """
    d = np.loadtxt(fname, delimiter=',')
    l = list()
    for i, row in enumerate(d):
        obsid = 'obshist' + str(i)
        lst = [obsid] + [mjd] + row.tolist()
        l.append(lst)
    return l
# main example :  create all the 'observations' separated by a day,

def obsMetaDataList(cadence=3.0, numepochs=20):
    """
    create a list of obsMetaData variables for a list of pointings.


    Parameters
    ----------
    cadence: float, optional, defaults to 3.0
        interval between epochs of observations
    numepochs: int, optional, defaults to 20
        number of epochs of observations


    Returns
    -------
    list of `~lsst.sims.catalogs.generation.db.ObservationMetaData` variables


    .. note: This will finally come from OpSims output, this is just a mock up. 
    """
    
    # List of MJDs with a three day cadence with 20 epochs 
    myMJDS = [570123.15 + cadence*i for i in range(numepochs)]
    filters = ['u', 'g', 'r', 'i', 'z', 'y']
    unrefRA = 5.0 
    unrefDec = 15.0 
    boundLen = 0.015
    boundType = 'circle'

    obsMetaDataList = []
    for mjd in myMJDS:
        obsMetaDataList.append(ObservationMetaData(boundType=boundType,
                               unrefractedRA=unrefRA,
                               unrefractedDec=unrefDec,
                               boundLength=boundLen,
                               bandpassName=filters,
                               mjd=mjd))
    return obsMetaDataList


def writeCatalogtoDB(dbfile, dbtable, ascii_root, galdb, obsMetaDataList):
    '''
    Write a set of instance catalogs to a sqlite3 database file called dbfile,
    deleting the file if it already exists. This is done in 2 steps: first the
    write method of instance catalogs is used to write to ASCII files that are
    not removed, and then these are read into the database. The ASCII files
    are not removed. It is assumed that the directory in which these files
    are written to exists.


    Parameters
    ----------
    dbfile : string
        absolute path to the database file
    dbtable : string
        table name in the database file
    ascii_root : str
        Determines the set of filenames that the instance catalogs are written
        to. The Filenames are of the form ascii_root_i.txt where i is integer
        (with no padding) the MJD in observations. For example, the instance
        catalog corresponding to the 3rd observation is written to a ASCII file
        of name 'data/SNIaCat_3.txt' relative to `pwd` if ascii_root is
        'data/SNIaCata_'
    galdb : `~CatalogDBObject.from_objid` galaxy catalog database  
    obsMetaDataList : list of
        `~lsst.sims.catalogs.generation.db.ObservationMetaData` variables


    Returns
    -------
    None


    ..note: But output is written to disk
    '''

    # erase database if it already exists
    if os.path.exists(dbfile):
        print 'deleting previous database'
        os.unlink(dbfile)

    # Setup connection to write
    connection = sqlite3.connect(dbfile)
    curs = connection.cursor()
    curs.execute('CREATE TABLE if not exists mysncat \
            (id TEXT, mjd FLOAT, snid INT, snra FLOAT, sndec FLOAT,\
            z FLOAT, t0 FLOAT, c FLOAT, x1 FLOAT, x0 FLOAT,\
            flux_u FLOAT, flux_g FLOAT, flux_r FLOAT, flux_i FLOAT,\
            flux_z FLOAT, flux_y FLOAT,\
            mag_u FLOAT, mag_g FLOAT, mag_r FLOAT, mag_i FLOAT,\
            mag_z FLOAT, mag_y FLOAT)')
    # Catalog, and range over which it is written
    # myMJDS = [570123.15 + 3.*i for i in range(20)]
    # for i, myMJD in enumerate(myMJDS):
    for i, obsmetadata in enumerate(obsMetaDataList):
        myObsMD = obsmetadata 
        catalog = SNIaCatalog(db_obj=galDB,
                              obs_metadata=myObsMD)
        print "====================================="
        print i, type(catalog.usedlsstbands()), catalog.obs_metadata.mjd
        print "====================================="
        # fname = "data/SNIaCat_" + str(i) + ".txt"
        fname = ascii_root + str(i) + ".txt"
        catalog.write_catalog(fname)
        l = _file2lst(fname, i, mjd=catalog.obs_metadata.mjd)
        recs = sq.array2dbrecords(l)
        exec_str = sq.insertfromdata(tablename=dbtable, records=recs,
                                     multiple=True)
        curs.executemany(exec_str, recs)
    connection.commit()
    connection.close()


def getLCsFromDB(dbfile, dbtable, lc_root):
    """
    Obtain light curves from a database table dbtable in a database file dbfile
    and write them to light curve files with names given by lc_rootsnid.dat.

    Parameters
    ----------
    dbfile : string
        absolute path to the database file

    dbtable : string
        table name in the database file

    lc_root : string
        string prepended to the the SNID and .dat to create the file name of
        the light curve. For example, lc_root = 'data/SN_' will write the light
        curve of SN with snid '1234' to 'data/SN_1234.dat'

    Returns
    -------

    None

    ..note: This function writes to disk at lc_rootsnid.dat for each SN.
    """
    connection = sqlite3.connect(dbfile)
    df = pd.read_sql('SELECT * FROM ' + dbtable, connection)
    grouped = df.groupby('snid')
    snids = grouped.groups.keys()
    lcs = [] 
    for snid in snids:
        fname = lc_root + str(snid) + '.dat'
        mylc = df.loc[grouped.groups[snid]]
        lcs.append((snid, mylc))
        mylc.to_csv(fname, na_rep="NaN", index=False)
    connection.close()
    return lcs

#def convert2SNcosmo(lc):
    """
    convert a lc into SNCosmo photometric data format returning an
    `~astropy.Table`

    """
    # metadata in lc is snid, snra, sndec, z, t0, x1, x0, c
#    RA = sn.

if __name__ == "__main__":


    galDB = CatalogDBObject.from_objid('galaxyTiled')
    obsmetaDataList = obsMetaDataList(cadence=3.0, numepochs=20)
    writeCatalogtoDB(dbfile='data/sncat.db',
                     dbtable='mysncat',
                     ascii_root='data/SNIaCat_',
                     galdb=galDB,
                     obsMetaDataList=obsmetaDataList)
    getLCsFromDB(dbfile='data/sncat.db',
                 dbtable='mysncat',
                 lc_root='data/LightCurves/SN_')
