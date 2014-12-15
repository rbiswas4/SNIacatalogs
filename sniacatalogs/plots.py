#/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
data = np.loadtxt("SNIaCat.txt", delimiter = ", ")
bins = np.arange( data[:,4].min(), data[:,4].max(), 100.)
plt.hist(data[:,4], histtype='step', bins=100.)
plt.axvline(570000.)
plt.show()
plt.hist(data[:,6], histtype="step")
plt.xlabel("c")
plt.show()
plt.hist(data[:,5], histtype="step")
plt.xlabel("x1")
plt.show()
