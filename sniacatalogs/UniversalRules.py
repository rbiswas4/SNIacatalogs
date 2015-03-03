#!/usr/bin/env python

import numpy as np


class Universe(object):
    """
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
        vals = np.zeros(shape=(self.numobjs, 4))

        for i, v in enumerate(vals):
            t0val = self.drawFromT0Dist() 
            
    def drawFromT0Dist(self, hostid):
        '''

        '''
        # Will not use hostid for now, but this is there so that 
        # later on one could obtain z, hostmass etc.
        t0val = np.uniform(-hundredyear / 2.0 + self.midSurveyTime, 
                           hundredyear / 2.0 + self.midSurveyTime)
        if self.suppressDimSN:
            if np.abs(t0val - self.obs_metadata.mjd) > self.maxTimeSNVisible:
                t0val = self.badvalues

