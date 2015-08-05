#!/usr/bin/env python

import numpy as np
import os
import sqlite3
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

import sniacatalogs as snia
from sniacatalogs import sncat
import lightcurve_utils as lu

from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catUtils.mixins import CosmologyMixin
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
import eups


print sncat.__file__



def samplePatchOnSphere(phi, theta, delta, size):
    '''
    (not) Uniformly distributes samples on a spherical patch between phi \pm delta,
    and theta \pm delta.
    
    will modify to fix this later
    Parameters
    ----------
    phi: float, mandatory, radians
        center of the spherical patch in ra with range 
    theta: float, mandatory, radians
    delta: float, mandatory, radians
    size: int, mandatory
        number of samples
    '''
    u = np.random.uniform(size=size)
    v = np.random.uniform(size=size)

    # phivals = delta * (2. * u - 1) + phi
    phivals = 2. * delta* u + (phi - delta )
    phivals = np.where ( phivals >= 0., phivals, phivals + 2. * np.pi)
    
    # use conventions in spherical coordinates
    theta = np.pi/2.0 - theta
    
    thetamax = theta + delta
    thetamin = theta - delta
    # CDF is cos(thetamin) - cos(theta) / cos(thetamin) - cos(thetamax)
    a = np.cos(thetamin) - np.cos(thetamax)
    thetavals = np.arccos(-v * a + np.cos(thetamin))
    # Get back to -pi/2 to pi/2 range of decs
    thetavals = np.pi/2.0 - thetavals 
    return phivals, thetavals


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
    phi = np.radians(mydict['unrefractedRA'])
    theta = np.radians(mydict['unrefractedDec'])
    equalrange = np.radians(mydict['boundLength'])
    ravals, thetavals = samplePatchOnSphere(phi=phi, theta=theta, delta=equalrange, size=size)
    return ravals, thetavals


def _createFakeGalaxyDB(dbname, ObsMetaData, size=10000, seed=1):
    '''
    Create a local sqlite galaxy database having filename dbname with variables
    id, raJ2000, decJ2000 and redshift, having number of rows =size, and having
    overlap with ObsMetaData.
    '''
    cleanDB(dbname)
    conn = sqlite3.connect(dbname)
    curs = conn.cursor()
    curs.execute('CREATE TABLE if not exists gals (id INT, raJ2000 FLOAT,                  decJ2000 FLOAT, redshift FLOAT)')

    np.random.seed(seed)
    samps = sample_obsmetadata(ObsMetaData, size=size)

    for count in range(size):
        id = 1000000 + count

        # Main Database should have values in degrees
        ra = np.degrees(samps[0][count])
        dec = np.degrees(samps[1][count])
        redshift = np.random.uniform()
        row = tuple([id, ra, dec, redshift])
        exec_str = insertfromdata(tablename='gals', records=row,
                                     multiple=False)
        curs.execute(exec_str, row)

    conn.commit()
    conn.close()
    return samps


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


class myGalaxyCatalog(CatalogDBObject):
    '''
    Create a like CatalogDBObject connecting to a local sqlite database
    '''

    objid = 'mytestgals'
    tableid = 'gals'
    idColKey = 'id'
    objectTypeId = 0
    appendint = 10000
    database = 'testdata/galcat.db'
    raColName = 'raJ2000'
    decColName = 'decJ2000'
    driver = 'sqlite'

    # columns required to convert the ra, dec values in degrees
    # to radians again
    columns = [('id', 'id', int),
               ('raJ2000','raJ2000 * PI()/ 180. '),
               ('decJ2000','decJ2000 * PI()/ 180.'),
               ('redshift', 'redshift')]


### Now create a fake galaxy database

# Using the functions above, we will create a fake galaxy database with galaxies in a small region of sky (tuned for the opsims areas that will be searched for later), spatially distributed roughly uniformly 


obsMetaDataforCat = ObservationMetaData(boundType='circle',
                                          boundLength=np.degrees(0.25),
                                          unrefractedRA=np.degrees(0.13),
                                          unrefractedDec=np.degrees(-1.2),
                                          bandpassName=
                                          ['u', 'g', 'r', 'i', 'z', 'y'],
                                          mjd=49350.)



sample_obsmetadata(obsMetaDataforCat, size=1)

obsMetaDataforCat.boundLength

obsMetaDataforCat._boundLength


# In[21]:

obsMetaDataforCat.bounds.radius


# In[22]:

obsMetaDataforCat.boundType


vals = _createFakeGalaxyDB(dbname='testData/galcat.db',
                    ObsMetaData=obsMetaDataforCat,
                    size=1000000,
                    seed=1)


# In[24]:

vals


# In[25]:

plt.plot(vals[0][:1000], vals[1][:1000], '.')
plt.axvline(2. * np.pi, color='k', lw=2.)
plt.axvline(0., color='k', lw=2.)
plt.axhline(np.pi, color='k', lw=2.)
plt.axhline(-np.pi, color='k', lw=2.)
plt.plot([0.13], [-1.2], 'rs', markersize=8)


# Fake galaxycatalog object 

# Galaxies (blue points) simulated in a square patch of RA, DEC around the red point

# In[28]:

galDB = myGalaxyCatalog()


# Check that there really are galaxies in the catalog by writing an instance catalog of the galaxies to disk

# In[29]:

class galCopy(InstanceCatalog):
    column_outputs = ['id', 'raJ2000', 'decJ2000', 'redshift']
    override_formats = {'raJ2000': '%8e', 'decJ2000': '%8e'}


# In[30]:

galphot = galCopy(db_obj=galDB, obs_metadata=obsMetaDataforCat)


# In[31]:

galphot.write_catalog('gals.dat')


# In[32]:



# The next step takes a set of mjds and pointings made up here (ie. this should really come from Opsims, but to simplify this notebook, and build a light curve, we will assume that we decided to point at the same location everyday. And to first try out
# an instance catalog we can read easily, and understand, we will first do one day.
# 
# 
# So, let us instantiate an ObservationMetaData object, and write the catalog associated with it

### Building Light Curves from Instance Catalogs

# We would like to create a number of supernova instance catalogs and then build the light curves from the catalogs. To do this correctly, we would like to use the `observation_metadata` associated with a number of conscutive OpSIM pointings. 

# In[34]:

opsimPath = os.path.join(eups.productDir('sims_data'),'OpSimData')
opsimDB = os.path.join(opsimPath,'opsimblitz1_1133_sqlite.db')

# from Tuscon AHM notebook from Scott
# This OPSIM DB is provided in sims_data. This creates a list of opsim pointings
# that I have checked. This is a tuned notebook
generator = ObservationMetaDataGenerator() #database = opsimPath, driver='sqlite')
obsMetaDataResults = generator.getObservationMetaData(limit=100,
                                                      fieldRA=(5.0, 8.0), 
                                                      fieldDec=(-85.,-60.),
                                                      expMJD=(49300., 49400.),
                                                      boundLength=0.015,
                                                      boundType='circle')


# In[35]:

# How many pointings do we have? 
print len(obsMetaDataResults)


# What are the RA, DEC values of the pointings? And how do they compare with the positions of galaxies we created in the database?

# In[36]:

def coords(x): 
    return np.radians(x.summary['unrefractedRA']), np.radians(x.summary['unrefractedDec'])


# In[37]:

v = zip(*map(coords, obsMetaDataResults))


# In[38]:

plt.plot(v[0], v[1], 'ko', markersize=4)
plt.plot(vals[0][:1000], vals[1][:1000], '.')
plt.axvline(2. * np.pi, color='k', lw=2.)
plt.axvline(0., color='k', lw=2.)
plt.axhline(np.pi, color='k', lw=2.)
plt.axhline(-np.pi, color='k', lw=2.)
plt.plot(v[0], v[1], 'ko', markersize=8)


# The few black points represent the opsim pointings, the blue points represent the locations of galaxies we created. This shows that there is overlap in space.

# We now want to write the light curves to a database. The way this function actually works is it creates the instance catalogs 
# for each of the obs_metadata points, writes them to disk using the `InstanceCatalog.write_catalog` method, reads the files from disk (right now the filenames are chosen in the same order as obs_metadata objects), and writes them to a sqlite database.

# In[39]:

if not os.path.exists('data/NewLightCurves'): 
    os.makedirs('data/NewLightCurves')


#### The Supernova Catalog

# In[40]:

column_outputs=['flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y', 'mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y']


# In[41]:

catalog = sncat.SNIaCatalog(db_obj=galDB, 
                            obs_metadata=obsMetaDataResults[0], 
                            column_outputs=['flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y', 'mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y'])


# In[42]:

catalog.write_catalog('test_SN_why.dat')


# In[42]:

#Set some 


# In[43]:

catalog.midSurveyTime = 49350.
catalog.averageRate = 1.


# In[44]:

catalog.write_catalog('tmp.dat')


# In[45]:



# In[46]:

f = catalog.get_SNsed()


# In[47]:

for sed in f:
    plt.plot(sed.wavelen, sed.flambda)


# In[48]:



# In[49]:

lu.writeCatalogtoDB(dbfile='data/sncat_real.db', dbtable='mysncat', ascii_root='data/SNIa_real_', galdb=galDB, 
                    obsMetaDataList=obsMetaDataResults, averageRate=1., midSurveyTime=49350., column_outputs=column_outputs)
    


# The following function then queries the database, and groups the data by a SNID number to build the light curves.

# In[50]:

lcs = lu.getLCsFromDB(dbfile='data/sncat_real.db',
                 dbtable='mysncat',
                 lc_root='data/LightCurves/SN_real_')


# We can check how many light curves are there, and we can display each light curve in the following way

# In[51]:

len(lcs)


# In[52]:

lcs[0][1]


# In[53]:

type(lcs[3][1])


# In this case, we have light curves associated with two different objects, and we can now plot them.

# In[54]:

fig, ax = plt.subplots(figsize=(24,8))
lu.plotlc(lcs[0])
lu.plotlc(lcs[1], marker='sr')
ax.invert_yaxis()


# The circles and squares are two different SN which were in lcs. The colors represent the different bands.

# In[55]:

epsilon = 0.01
theta = -np.pi/4
np.cos(theta) - np.cos(theta - epsilon)


# In[55]:

import os
from lsst.sims.GalSimInterface import GalSimGalaxies, GalSimSN
from lsst.sims.photUtils import PhotometricParameters
# from HackTogetherNIRCam import ReturnNIRCam


# In[56]:

class GalsimSNCatalog(GalSimSN):

    photParams = PhotometricParameters()
    
    



import galsim
print galsim.__file__
from lsst.sims.GalSimInterface import PSFbase

class myCrazyPSF(PSFbase):
    wavelength_dependent = False
    
    def _getPSF(self, xPupil=None, yPupil=None, **kwargs):
        gaussian = galsim.Gaussian(sigma=1.0)
        xp = np.abs(xPupil)
        yp = np.abs(yPupil)
        if xp is None:
            xp is 0.0
        if yp is None:
            yp = 0.0
            
        if xp<yp:
            minor = xp
            major = yp
        else:
            minor = yp
            major = xp

        _psf = gaussian.shear(q=(minor+1.0)/(major+1.0), beta=(xPupil+yPupil)*galsim.radians)
        return _psf


# In[58]:

class CrazyPSFCatalog2(GalsimSNCatalog):
    PSF = myCrazyPSF()


# In[59]:

cat = CrazyPSFCatalog2(db_obj=galDB,  obs_metadata=obsMetaDataResults[0])

cat.midSurveyTime = 49350.
cat.averageRate = 1.0
print 'Starting write'
cat.write_catalog('galsim.dat')
print 'Starting write images'
cat.write_images('test.fits')




