import sqlite3
import os
import numpy as np
import matplotlib.pyplot as plt
import sncosmo
from astropy.units import Unit

def getlsstbandpassobjs(loadsncosmo=True, loadcatsim=True, plot=False):
    """
    General utility to return a list of the baseline LSST bandpasses loaded
    as catsim bandpass objects, and register them as SNCosmo bandpasses
    accessible through strings like 'LSSTu'.
    args:
        loadsncosmo: Bool, optional, defaults to True
            variable to decide whether to register the LSST bandpasses as
            SNCosmo registered bandpass objects accessible through strings
            like 'LSSTu'
        loadcatsim : Bool, optional, defaults to True
            variable to decide whether to set up catsim bandpass objects
            are return the list of u,g,r,i,z,y bandpasses
    returns:
        if loadcatsim is true, list of catsim bandpass objects corresponding
        to LSST baseline u, g, r, i, z, y filters.
        if loadcatsim is False, return is None

    Examples:

    """
    bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
    banddir = os.path.join(os.getenv('THROUGHPUTS_DIR'), 'baseline')
    lsstbands = []
    lsstbp = {}

    for band in bandPassList:
        # setup sncosmo bandpasses
        bandfname = banddir + "/total_" + band + '.dat'
        if loadsncosmo:
            # register the LSST bands to the SNCosmo registry
            # Not needed for LSST, but useful to compare independent codes
            # Usually the next two lines can be merged,
            # but there is an astropy bug currently which affects only OSX.
            numpyband = np.loadtxt(bandfname)
            sncosmoband = sncosmo.Bandpass(wave=numpyband[:, 0],
                                           trans=numpyband[:, 1],
                                           wave_unit=Unit('nm'),
                                           name='LSST' + band)
            sncosmo.registry.register(sncosmoband, force=True)
        if loadcatsim:
            # Now load LSST bandpasses for catsim
            lsstbp[band] = Bandpass()
            lsstbp[band].readThroughput(bandfname, wavelen_step=wavelenstep)
            lsstbands.append(lsstbp[band])
    fifilterfigs, filterax = plt.subplots()
    if plot:
        for band in bandPassList:
            b = sncosmo.get_bandpass('LSST' + band)
            filterax.plot(b.wave, b.trans, '-k', lw=2.0)
            filterax.set_xlabel(r'$\lambda$, ($\AA$)')
            filterax.set_ylabel(r'transmission')
            filterax.set_label(r'LSST Filter Transmission Functions')
        plt.show()

    if loadcatsim:
        return lsstbands
    else:
        return None

#
#                             SQLLITE UTILS
#

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

def samplePatchOnSphere(phi, theta, delta, size):
    '''
    Uniformly distributes samples on a spherical patch between phi \pm delta,
    and theta \pm delta.
    Parameters
    ----------
    phi: float, mandatory
        center of the spherical patch
    theta: float, mandatory

    delta: float, mandatory

    size: int, mandatory
        number of samples

    '''
    u = np.random.uniform(size=size)
    v = np.random.uniform(size=size)

    phivals = delta * (2. * u - 1) + phi

    # thetavals
    thetamax = theta + delta
    thetamin = theta - delta
    a = np.cos(thetamax) - np.cos(thetamin)
    thetavals = np.arccos(v* a + np.cos(thetamin))
    return phivals, thetavals


def sample_obsmetadata(obsmetadata, size=1):
    '''
    Sample a square patch on the sphere overlapping obsmetadata
    field of view by picking the area enclosed in
    obsmetadata.unrefractedRA \pm obsmetadata.boundLength
    obsmetadata.unrefractedDec \pm obsmetadata.boundLength

    Parameters
    ----------
    obsmetadata: instance of
        `sims.catalogs.generation.db.ObservationMetaData`

    size: integer, optional, defaults to 1
        number of samples


    Returns
    -------

    tuple of ravals, decvalues
    '''
    mydict = obsmetadata.summary
    phi = mydict['unrefractedRA']
    theta = mydict['unrefractedDec']
    equalrange = mydict['boundLength']
    ravals, thetavals = samplePatchOnSphere(phi=phi, theta=theta, delta=equalrange, size=size)
    return ravals, thetavals


def cleanDB(dbname, verbose=True):
    '''
    Deletes the database dbname from the disk.
    Parameters
    ----------
    dbname: string, mandatory
        name (abs path) of the database to be deleted
    verbose: Bool, optional, defaults to True

    '''

    if os.path.exists(dbname):
        if verbose:
            print "deleting database ", dbname
        os.unlink(dbname)
    else:
        if verbose:
            print 'database ', dbname, ' does not exist'


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
