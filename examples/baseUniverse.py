#!/usr/bin/env python

from UniversalRules import Universe
import numpy as np

class SNUniverse( Universe):

    def __init__(self):

        self.numobjs = 1
        self.suppressHighzSN = True
        self.badvalues = np.nan
        self.midSurveyTime = 570000.
        self.suppressDimSN = False
        self.mjdobs = 571203.15
        self.maxTimeSNVisible = 100.

        return 

if __name__ == '__main__':

    # from UniversalRules import Universe
    snu = SNUniverse()

    # class myU(erse, Universe):
    #    pass

    # uu = myU(snu)

    print snu.numobjs
    print snu.suppressDimSN 
    print snu.suppressHighzSN
    print snu.badvalues
    print snu.mjdobs
    print snu.maxTimeSNVisible
    np.random.seed(144)
    x1 = snu.drawFromx1Dist()
    c = snu.drawFromcDist()
    t0 = snu.drawFromcDist()
    print snu.drawFromcDist()
    print 'PARAMDIST', snu.SNparamDistfromHost(hostz=[1.0], hostid=[4], hostmu=[42.])
