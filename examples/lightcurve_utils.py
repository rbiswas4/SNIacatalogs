#!/usr/bin/env python

import numpy as np
import os
import sqlite3
import pandas as pd
from sniacatalogs import sncat
import matplotlib.pyplot as plt

def _file2lst(fname, i, mjd):
    """
    creates a lst of observations from a Catalog written out
    """
    d = np.loadtxt(fname, delimiter=',', ndmin=2)
    l = list()
    if len(d) == 0:
        return l
    for i, row in enumerate(d):
        obsid = 'obshist' + str(i)
        lst = [obsid] + [mjd] + row.tolist()
        l.append(lst)
    return l


def array2dbrecords(data):
    """
    takes an iterable such as an array, structured array or list and converts it into a a list of tuples suitable
    for insertion into a sqlite 3 database. 
    args:
        data: mandatory,
            list or input iterable
    returns:
        suitable list of tuples
        
    """
    l = list()
    for row in data:
        l.append(tuple(row))
    return l

def insertfromdata(tablename, records, multiple=True):
    """
    construct string to insert multiple records into sqlite3 database 
    args:
        tablename: str, mandatory
            Name of table in the database. 
        records: set of records
        multiple:
    returns:
    """
    if multiple:
        lst = records[0]
    else:
        lst = records
    s = 'INSERT INTO ' + str(tablename) + ' VALUES '
    s += "( " + ", ".join(["?"]*len(lst)) + ")"
    return s

def obsMetaDataList(startdate=570180, cadence=3.0, numepochs=20):
    """
    create a list of obsMetaData variables for a list of pointings.


    Parameters
    ----------
    startdate: float, optional, defaults to 570180
        starting date of sequence
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
    myMJDS = [startdate + cadence*i for i in range(numepochs)]
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


def writeCatalogtoDB(dbfile, dbtable, ascii_root, galdb, obsMetaDataList,
                     midSurveyTime=5700., averageRate=36500.):
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
        absolute path to the database file to write to
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
    midSurveyTime: float, optional defaults to 57000.
        Set the midSurveyTime of the SN Catalog
    averageRate: float, optional defaults to 36500.
        The SN goes of with a probability of 1/averageRate every day in a galaxy


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
    curs.execute('CREATE TABLE if not exists mysncat             (id TEXT, mjd FLOAT, snid INT, snra FLOAT, sndec FLOAT,            z FLOAT, t0 FLOAT, c FLOAT, x1 FLOAT, x0 FLOAT,            flux_u FLOAT, flux_g FLOAT, flux_r FLOAT, flux_i FLOAT,            flux_z FLOAT, flux_y FLOAT,            mag_u FLOAT, mag_g FLOAT, mag_r FLOAT, mag_i FLOAT,            mag_z FLOAT, mag_y FLOAT)')
    # Catalog, and range over which it is written
    # myMJDS = [570123.15 + 3.*i for i in range(20)]
    # for i, myMJD in enumerate(myMJDS):
    for i, obsmetadata in enumerate(obsMetaDataList):
        myObsMD = obsmetadata
        catalog = sncat.SNIaCatalog(db_obj=galdb,
                              obs_metadata=myObsMD)
        catalog.averageRate = averageRate
        catalog.midSurveyTime = midSurveyTime
        print "====================================="
        print i,  catalog.obs_metadata.mjd
        print "====================================="
        # fname = "data/SNIaCat_" + str(i) + ".txt"
        fname = ascii_root + str(i) + ".txt"
        print catalog.averageRate, catalog.midSurveyTime
        catalog.write_catalog(fname)
        l = _file2lst(fname, i, mjd=catalog.obs_metadata.mjd)
        # print l
        if len(l) > 0:
            recs = array2dbrecords(l)
            exec_str = insertfromdata(tablename=dbtable, records=recs,
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



def plotlc(lc,marker='o'):
    #fig, ax = plt.subplots()
    ax = plt.gca()
    ax.plot(lc[1]['mjd'], lc[1]['mag_u'], marker, color='b')
    ax.plot(lc[1]['mjd'], lc[1]['mag_g'], marker, color='g')
    ax.plot(lc[1]['mjd'], lc[1]['mag_r'], marker, color='r')
    ax.plot(lc[1]['mjd'], lc[1]['mag_i'], marker, color='m')
    ax.plot(lc[1]['mjd'], lc[1]['mag_z'], marker, color='y')
    ax.plot(lc[1]['mjd'], lc[1]['mag_y'], marker, color='k')
    #ax.invert_yaxis()
    ax.set_xlabel('MJD')
    ax.set_ylabel('mag')
    



