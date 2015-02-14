"""
Class describing the SN object itself. The SN object derives from
SNCosmo.Model and provides a model sed of SNIa based on the SALT2 model-like
 model often called 'salt2-extended'. This model is extended to longer 
 wavelength ranges compared to models in Guy10 or Betoule14, at the cost of
 larger model varianc. SNObject has additional attributes such as ra, dec. and
additional methods to calculate band magnitudes using the LSST software stack
after applying MW extinction:
 -  calc_mags which use the magnitude calculations in LSST stack
 -  extinction which use the extinction calculations in LSST stack

"""
import numpy as np
import os

from lsst.sims.photUtils.Photometry import PhotometryStars, Sed, Bandpass
from lsst.sims.photUtils.EBV import EBVbase

import sncosmo


wavelenstep = 0.1
class SNObject (sncosmo.Model):
    """
    Extension of the SNCosmo `TimeSeriesModel` to include more parameters and
    use methods in the catsim stack. We constrain ourselves to the use of a
    specific SALT model for the Supernova (Salt2-Extended), and set its MW
    extinction to be 0, since we will use the LSST software to calculate
    extinction.

    Parameters
    ----------
    ra : float
         ra of the SN in degrees
    dec : float
        dec of the SN in degrees

    
    Attributes
    ----------
    ra : float or None 
        ra of the SN in radians
        
    dec : float or None
        dec of the SN in radians

    skycoord: `np.ndarray' of size 2 or None
        np.array([[ra], [dec]]), which are in radians

    ebvofMW: float or None
        mwebv value calculated from the self.skycoord if not None, or set to 
        a value using self.set_MWebv. If neither of these are done, this value
        will be None, leading to exceptions in extinction calculation.
        Therefore, the value must be set explicitly to 0. to get unextincted
        quantities.

        
    Methods
    -------
    """
    def __init__(self, ra=None, dec=None):
        """
        Instantiate class

        Parameters
        ----------

        ra : float
            ra of the SN in degrees
        dec : float
            dec of the SN in degrees

        """
        dust = sncosmo.CCM89Dust()
        sncosmo.Model.__init__(self, source="salt2-extended",
                       effects=[dust, dust], effect_names=['host', 'mw'],
                       effect_frames=['rest', 'obs'])
        # Current implementation of Model has a default value of mwebv = 0.
        # ie. no extinction, but this is not part of the API, so should not
        # depend on it, set explicitly in order to unextincted SED from
        # SNCosmo. We will use catsim extinction from photUtils.
        self.set(mwebv=0.)
        # If we know ra, dec in degree
        self.ra = ra
        self.dec = dec
        self.skycoord = None
        if self.dec is not None:
            self.dec = dec * np.pi / 180.0 
        if self.ra is not None:
            self.ra = ra * np.pi / 180.0

        self.lsstmwebv = EBVbase()
        self.ebvofMW = None
        self.mwEBVfromMaps()
        self._seed = None
        return

    # Comment to self: Don't see why having the seed is of any help here
    # will probably remove
    # @property
    # def seed(self):
    #    return self._seed

    def set_MWebv(self, value):
        """
        if mwebv value is known, this can be used to set the attribute
        _mwebv of the SNObject class to the value (float).

        Parameters
        -----------

        value : float, mandatory
                value of mw extinction parameter E(B-V) in mags to be used in
                applying extinction to the SNObject spectrum

        Returns
        -------
        None

        .. note:: For a large set of SN, one may use fast `np.ndarray` valued 
                  functions to obtain an array of such values, and then set 
                  the values from such an array.
        """
        #self._mwebv = value
        self.ebvofMW = value
        return

    def mwEBVfromMaps(self):
        """
        set the attribute _mwebv of the class from the ra and dec
        of the SN. If the ra or dec attribute of the class is None,
        set this attribute to None.

        Parameters
        ----------
        None

        Returns 
        -------
        None

        .. note:: This function must be run after the class has attributes ra 
                  and dec set. In case it is run before this, the mwebv value 
                  will be set to None.

        """

        # If ra or dec is None, leave mwebv as initialized 
        if self.ra is None or self.dec is None:
            return
        # else set skycoord
        self.skycoord = np.array([[self.ra], [self.dec]])
        self.ebvofMW  = self.lsstmwebv.calculateEbv(equatorialCoordinates=
                                                  self.skycoord)[0]
        return

    def bandMags(self, bandpassobjects, time, phiarray=None):
        """
        return a numpy array of magnitudes of the SN spectrum in the ab
        magnitude system.

        Parameters
        ----------
        bandpassobjects: mandatory, list of bandpass objects
                         a list of LSST catsim bandpass objects
        time: mandatory, float
              MJD at which this is evaluated

        Returns
        -------
        `np.ndarray` of mag values for each band in lsstbandpass. 

        .. note:: Unphysical values of the flux density are reported as `np.nan`
        """

        # if phiarray is None:
        if self.parameters[0] > 1.2:
            return [np.nan]*len(bandpassobjects)

        filterwav = bandpassobjects[0].wavelen
        # else: 
        #    filterwav = phiarray[0].wavelen
        # print 'new maccalc'
        # print filterwav.max(), filterwav.min()
        # print self.minwave(), self.maxwave(), self.parameters
        
        SEDfromSNcosmo = Sed(wavelen=filterwav,
                             flambda=self.flux(time=time,
                                               wave=filterwav*10.)*10.)
        # Check if I need this sync
        # SEDfromSNcosmo.synchronizeSED(wavelen_min=filterwav[0],
        #                              wavelen_max=filterwav[-2],
        #                              wavelen_step=wavelenstep)
        # Apply LSST extinction
        ax, bx = SEDfromSNcosmo.setupCCMab()
        SEDfromSNcosmo.addCCMDust(a_x=ax, b_x=bx, ebv=self.ebvofMW)
        if phiarray is None:
            phiarray, dlambda = SEDfromSNcosmo.setupPhiArray(bandpassobjects)
        SEDfromSNcosmo.synchronizeSED(wavelen_min=filterwav[0],
                                      wavelen_max=filterwav[-2],
                                      wavelen_step=wavelenstep)
        
        return SEDfromSNcosmo.manyMagCalc(phiarray, wavelen_step=wavelenstep)
