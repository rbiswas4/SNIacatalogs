#!/usr/bin/env python

import numpy as np
import os
import timeit

from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import ObservationMetaData
import lsst.sims.catUtils.baseCatalogModels as bcm


fname = 'results.dat'
galDB = CatalogDBObject.from_objid('galaxyTiled')

# Instance catalog class
class galCopy(InstanceCatalog):
    column_outputs = ['id', 'raJ2000', 'decJ2000', 'redshift']
    override_formats = {'raJ2000': '%8e', 'decJ2000': '%8e'}


# In[5]:

myMJD = 571203.15
myObsMD = ObservationMetaData(boundType='circle',
                                  boundLength=0.015,
                                  unrefractedRA=5.0,
                                  unrefractedDec=15.0,
                                  site='LSST',
                                  bandpassName=['u', 'g', 'r', 'i', 'z', 'y'],
                                  mjd=myMJD)


blarray = 1.75/200. * np.array(range(100, 199))
blarray = blarray [::5]
# blarray = np.array([0.0015, 0.015, 0.02])
# blarray

results = []
if os.path.exists(fname):
    os.remove(fname)
for bl in blarray:
    myObsMD = ObservationMetaData(boundType='circle',
                                  boundLength=bl,
                                  unrefractedRA=5.0,
                                  unrefractedDec=15.0,
                                  site='LSST',
                                  bandpassName=['u', 'g', 'r', 'i', 'z', 'y'],
                                  mjd=myMJD)
    gals = galCopy(db_obj=galDB, obs_metadata=myObsMD)
    timer = timeit.Timer('gals.write_catalog("gals.dat")', setup='from __main__ import gals' )
    res  = timer.repeat(repeat=3, number=1 )
    lines = sum(1 for _ in open('gals.dat'))
    it = [bl] + res + [lines]
    print it
    with open('results.dat', 'a') as f:
        strwrite = ' '.join(map(str, it)) + '\n'
        f.write(strwrite)
       

    results.append(it)
    print bl,  res[0]
    



# get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt
# import seaborn as sns
# sns.set()



results = np.array(results)


def factor(patchlen, maxlen=1.75  ):
    return (maxlen / patchlen)**2



fig, ax = plt.subplots(2,2)
ax[0, 0].errorbar(results[:, 0], np.mean(results[:, 1:4], axis=1), np.std(results[:, 1:4], axis=1), fmt='o')
ax[1, 0].errorbar(results[:, 0], np.mean(results[:, 1:4], axis=1)*factor(results[:, 0]), np.std(results[:, 1:4], axis=1)*factor(results[:, 0]), fmt='o')
ax[0, 1].errorbar(np.log10(results[:, -1]), np.mean(results[:, 1:4], axis=1), np.std(results[:, 1:4], axis=1), fmt='o')
ax[1, 1].errorbar(np.log10(results[:, -1]), np.mean(results[:, 1:4], axis=1)*factor(results[:, 0]), np.std(results[:, 1:4], axis=1)*factor(results[:, 0]), fmt='o')

import sys

fig.savefig('highend.pdf')
np.savetxt('highend.txt', results)
# In[35]:




# In[ ]:



