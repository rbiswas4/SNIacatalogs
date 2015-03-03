#!/usr/bin/env python

import numpy as np


class Universe(object):

    def __init__(self) :
        pass

    def SNCoordinatesFromHost(hostra, hostdec, hostz, numhosts):
        '''
        Distribution of SN coordinates and velocities given a set of host
        coordinates and velocities.
    
        '''
        _sndec = hostdec
        _snra = hostra
        snz = hostz
        snvra = np.zeros(numhosts)
        snvdec = np.zeros(numhosts)
        snvr = np.zeros(numhosts)
    
        return _snra, _sndec, snz , snvra, snvdec, snvr 
