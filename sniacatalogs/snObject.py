"""
Class describing the SN object itself. The SN object derives from
SNCosmo.Model and provides a sed of SNIa based on the SALT2 model-like
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

    Examples
    --------
    >>> SNObject  = SNObject(ra=30., dec=60.)
    >>> SNObject._ra
    >>> 0.5235987755982988
    >>> SNObject._dec
    >>> 1.0471975511965976
    """
    def __init__(self, ra=None, dec=None):
        """
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

        # ra and dec passed as parameters are in degrees
        self._ra = ra
        self._dec = dec

        # ra and dec if passed as non Nne variables are converted into
        # radians
        if self._dec is not None:
            self._dec = dec * np.pi / 180.0
        if self._ra is not None:
            self._ra = ra * np.pi / 180.0

        self.lsstmwebv = EBVbase()
        self.ebvofMW = None
        if self._ra is not None and self._dec is not None:
            self.mwEBVfromMaps()
        return
    def summary(self):
        '''
        summarizes the current state of the SNObject class in a returned
        string.

        Parameters
        ----------
        None

        Returns
        -------
        Summary State in string format

        Examples
        --------
        >>> t = SNObject()
        >>> print t.summary()


        '''
        state = '  SNObject Summary      \n'

        state += 'Model = ' + '\n'
        state += 'z = '+ str(self.get('z')) + '\n'
        state += 'c = '+ str(self.get('c')) + '\n'
        state += 'x1 = '+ str(self.get('x1')) + '\n'
        state += 'x0 = '+ str(self.get('x0')) + '\n'
        state += 't0 = '+ str(self.get('t0')) + '\n'
        state += 'ra = '+ str(self._ra) + ' in radians \n'
        state += 'dec = '+ str(self._dec) + ' in radians \n'
        state += 'MW E(B-V) = ' + str(self.ebvofMW) + '\n'
        state += '+++++++++++++++++++++++\n'

        return state
    def setCoords(self, ra, dec):
        """
        set the ra and dec coordinate of SNObject to values in radians
        corresponding to the given values in degrees

        Parameters
        ----------
        ra: float, mandatory
            the ra in degrees
        dec: float, mandatory
            dec in degrees

        Returns
        -------
        None

        Examples
        --------
        >>> t = SNObject()
        >>> t.setCoords(ra=30., dec=90.)
        >>> t._ra
        >>> 0.5235987755982988
        >>> t._dec
        >>> 1.0471975511965976
        """

        self._ra = ra * np.pi / 180.
        self._dec = dec * np.pi / 180.

        return
    def set_MWebv(self, value):
        """
        if mwebv value is known, this can be used to set the attribute
        ebvofMW of the SNObject class to the value (float).

        Parameters
        ----------
        value: float, mandatory
               value of mw extinction parameter E(B-V) in mags to be used in
               applying extinction to the SNObject spectrum

        Returns
        -------
        None

        Examples
        --------
        >>> t = SNObject()
        >>> t.set_MWebv(0.)
        >>> 0.

        .. note:: For a large set of SN, one may use fast `np.ndarray` valued
                  functions to obtain an array of such values, and then set
                  the values from such an array.
        """
        self.ebvofMW = value
        return

    def mwEBVfromMaps(self):
        """
        set the attribute ebvofMW of the class from the ra and dec
        of the SN. If the ra or dec attribute of the class is None,
        set this attribute to None.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Examples
        --------
        >>> t = SNObject()
        >>> t.setCoords(ra=30., dec=60.)
        >>> t.mwEBVfromMaps()
        >>> t.ebvofMW
        >>> 0.977767825127

        .. note:: This function must be run after the class has attributes ra
                  and dec set. In case it is run before this, the mwebv value
                  will be set to None.

        """

        if self._ra is None or self._dec is None:
            raise ValueError('Cannot Calculate EBV from dust maps if ra or dec is `None`')
            return
        # else set skycoord
        self.skycoord = np.array([[self._ra], [self._dec]])
        self.ebvofMW  = self.lsstmwebv.calculateEbv(equatorialCoordinates=
                                                    self.skycoord)[0]
        return

    def SNObjectSED(self, time, wavelen=None, bandpassobject=None,
                    applyExitinction=True):
        '''
        return a `sims_photutils.sed`  object from the SN model at the
        requested time and wavelengths (or wavelengths of a bandpassobject)
        with extinction from MW according to the SED extinction methods. If
        the sed is requested at times outside the validity range of the
        model, the flux density is returned as 0. If the time is within the
        range of validity of the model, but the wavelength range requested
        is outside the range, then the returned fluxes are np.nan outside
        the range, and the model fluxes inside

        Parameters
        ----------
        time: float
            time of observation
        wavelen: `np.ndarray` of floats, optional, defaults to None
            array containing wavelengths in nm
        bandpassobject: `sims_photUtils.bandpass` object or dict thereof,
            optional, defaults to `None`. Using the dict assumes that the
            wavelength sampling and range is the same for all elements of
            the dict.
            if provided, overrides wavelen input and the SED is
            obtained at the wavelength values native to bandpass
            object.


        Returns
        -------
        `sims_photutils.sed` object containing the wavelengths and SED
        values from the SN at time time in units of ergs/cm^2/sec/nm


        .. note: If both wavelen and bandpassobject are `None` then exception,
                 will be raised.
        Examples
        --------
        >>> sed = SN.SNObjectSED(time=0., wavelen=wavenm)
        >>> # Usual Sed methods and attributes can be accessed
        '''

        if wavelen is None and bandpassobject is None:
            raise ValueError('A non None input to either wavelen or\
                              bandpassobject must be provided')

        #if bandpassobject present, it overrides wavelen
        if bandpassobject is not None:
            if isinstance(bandpassobject, dict):
                firstfilter = bandpassobject.keys()[0]
                bp = bandpassobject[firstfilter]
            else:
                bp = bandpassobject
            # remember this is in nm
            wavelen =  bp.wavelen

        flambda = np.zeros(len(wavelen))

        # self.mintime() and self.maxtime() are properties describing
        # the ranges of SNCosmo.Model in time.

        # Set SED to 0 beyond the model phase range, will change this if
        # SNCosmo includes a more sensible decay later.
        if time > self.mintime() and time < self.maxtime():


            # If SNCosmo is requested a SED value beyond the model range
            # it will crash. Try to prevent that by returning np.nan for
            # such wavelengths. This will still not help band flux calculations
            # but helps us get past this stage.

            flambda = flambda * np.nan

            # Convert to Ang
            wave = wavelen * 10.
            mask1 = wave > self.minwave()
            mask2 = wave < self.maxwave()
            mask = mask1 & mask2
            wave = wave[mask]
            # flux density dE/dlambda returned from SNCosmo in
            # ergs/cm^2/sec/Ang, convert to ergs/cm^2/sec/nm
            flambda[mask] = self.flux(time=time, wave=wave)
            flambda[mask] = flambda * 10.

        SEDfromSNcosmo = Sed(wavelen=wavelen, flambda=flambda)

        if not applyExitinction:
            return SEDfromSNcosmo

        # Apply LSST extinction
        ax, bx = SEDfromSNcosmo.setupCCMab()
        if self.ebvofMW is None:
            # self.mwEBVfromMaps()
            # if self.ebvofMW is None:
            raise ValueError('ebvofMW attribute cannot be None Type and must be\
                             set by hand using set_MWebv before this stage, or\
                             by using setcoords followed by mwEBVfromMaps\n')

        SEDfromSNcosmo.addCCMDust(a_x=ax, b_x=bx, ebv=self.ebvofMW)
        return SEDfromSNcosmo


    def catsimBandMags(self, bandpassobject, time, phiarray=None,
                       observedBandPassInds=None):
        """
        return a numpy array of AB magnitudes in the bandpasses

        Parameters
        ----------
        bandpassobject: mandatory, bandpass objects or ordered dict thereof
            LSST Catsim bandpass object or ordered dict of bandpass objects
            in which the fluxes are calculated

        time: mandatory, float
              MJD at which this is evaluated

        phiarray: Optional, `sims.photometry.phiarray` object, defaults to None
            if present, is used in speeding up the flux calculation, if not,
            phiarray is calculated in this routine.

        observedBandPassInds: Optional, list of ints, defaults to None
            list of indices corresponding to the obserevd filters in the ordered
            dict bandpassobject. If present, fluxes are only returned for these
            bands.
        Returns
        -------
        `np.ndarray` of mag values for each band in lsstbandpass.

        Examples
        --------
        >>> t = SNObject()
        """


        SEDfromSNcosmo = self.SNObjectSED(time=time,
                                          bandpassobject=bandpassobject)

        if phiarray is None:
            phiarray, dlambda = SEDfromSNcosmo.setupPhiArray(bandpassobject)

        if isinstance(bandpassobject, dict):
            firstfilter=bandpassobject.keys()[0]
            wavelenstep = bandpassobject[firstfilter].wavelen_step
        else:
            wavelenstep = bandpassobject.wavelen_step
        SEDfromSNcosmo.flambdaTofnu()
        return SEDfromSNcosmo.manyMagCalc(phiarray,
                                          wavelen_step=wavelenstep, 
                                          observedBandPassInd=observedBandPassInds)

    def catsimBandFluxes(self, bandpassobject, time, phiarray=None,
                         observedBandPassInds=None):
        """
        return a numpy array of fluxes of the SN spectrum in the ab
        magnitude system

        Parameters
        ----------
        bandpassobject: mandatory, LSST catsim bandpass object or dict thereof
            LSST Catsim bandpass object or ordered dict of bandpass objects
            in which the fluxes are calculated

        time: mandatory, float
              MJD at which this is evaluated

        phiarray: Optional, `sims.photometry.phiarray` object, defaults to None
            if present, is used in speeding up the flux calculation, if not,
            phiarray is calculated in this routine.

        observedBandPassInds: Optional, list of ints, defaults to None
            list of indices corresponding to the obserevd filters in the ordered
            dict bandpassobject. If present, fluxes are only returned for these
            bands.
        Returns
        -------
        `np.ndarray` of mag values for each observed band in lsstbandpass.

        Examples
        --------
        >>> from lsst.sims.photUtils.Photometry import PhotometryBase
        >>> pbase = PhotometryBase()
        >>> pbase.loadBandpassesFromFiles()
        >>> pbase.setupPhiArray_dict() 

        .. note:: Unphysical values of the flux density are reported as
        `np.nan`
        """
        SEDfromSNcosmo = self.SNObjectSED(time=time,
                                          bandpassobject=bandpassobject)

        if phiarray is None:
            phiarray, dlambda = SEDfromSNcosmo.setupPhiArray(bandpassobject)


        if isinstance(bandpassobject, dict):
            firstfilter=bandpassobject.keys()[0]
            wavelenstep = bandpassobject[firstfilter].wavelen_step
        else:
            wavelenstep = bandpassobject.wavelen_step
        SEDfromSNcosmo.flambdaTofnu()

        return SEDfromSNcosmo.manyFluxCalc(phiarray, wavelen_step=wavelenstep,
                                           observedBandPassInd=observedBandPassInds)
