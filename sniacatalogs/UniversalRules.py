#!/usr/bin/env python

import numpy as np


class Universe(object):
    """
    The base class must have the following attributes:

    numobjs
    suppressHighzSN
    badvalues
    midSurveyTime
    suppressDimSN
    mjdobs
    maxTimeSNVisible

    Mixin Class for sncat providing methods
    """


    def SNCoordinatesFromHost(self, hostra, hostdec, hostz):
        '''
        Distribution of SN coordinates and velocities given a set of host
        coordinates and velocities.
        '''
        suppressHighzSN = self.suppressHighzSN 
        numhosts = self.numobjs 

        _sndec = hostdec
        _snra = hostra
        snz = hostz
        snvra = np.zeros(numhosts)
        snvdec = np.zeros(numhosts)
        snvr = np.zeros(numhosts)
    
        return _snra, _sndec, snz , snvra, snvdec, snvr 

    def SNparamDistfromHost(self, hostz, hostid, hostmu):
        '''
        Distribution of SN model parameters given their hosts
        '''
        print 'Got here'
        vals = np.zeros(shape=(self.numobjs, 4))

        for i, v in enumerate(vals):
            np.random.seed(hostid[i])
            t0val = self.drawFromT0Dist() 
            vals[i, 3] = t0val
            if t0val is self.badvalues:
                continue
            cval = self.drawFromcDist()
            x1val = self.drawFromx1Dist()
            x0val = self.drawFromX0Dist(x1val, cval, hostmu=hostmu[i])

            print cval, x1val, t0val
            vals[i, 0] = cval
            vals[i, 1] = x1val
            vals[i, 2] = x0val 
        
        return vals
                    
            
    def drawFromx1Dist(self, **hostParams):
        """
        """
        return np.random.normal(0. , 1.0) 

    def drawFromcDist(self, **hostParams):
        """
        """
        return np.random.normal(0. , 0.1) 
    
    def drawFromX0Dist(self, x1val, cval, hostmu, **hostParams):
        """
        """
        import snObject

        # First draw an absolute BessellB magnitude for SN
        mabs =  np.random.normal(-19.3, 0.3)
        mag = mabs + hostmu

        sn = snObject.SNObject()
        sn.set(x1=x1val, c=cval)
        sn.source.set_peakmag(mag, band='bessellb', magsys='ab')
        x0val = sn.get('x0')
        
        return x0val

    def drawFromT0Dist(self, **hostParams):
        '''
        Distribution function of the time of peak of SN


        '''

        # Will not use hostid for now, but this is there so that 
        # later on one could obtain z, hostmass etc. This is useful to obtain
        # z and host dependent SN rates
        hundredyear = 100. * 365.
        t0val = np.random.uniform(-hundredyear / 2.0 + self.midSurveyTime, 
                           hundredyear / 2.0 + self.midSurveyTime)

        if self.suppressDimSN:

            # SN far away from peak will have zero flux, but writing zero flux
            # objects to catalogs will mean each instance catalog will have
            # too many objects. So, this variable will always stay True (except
            # for testing purposes), while the variable that will be changed
            # is maxTimeSNVisible

            if np.abs(t0val - self.mjdobs) > self.maxTimeSNVisible:
                t0val = self.badvalues
        return t0val
