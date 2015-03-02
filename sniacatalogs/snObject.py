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

        # ra and dec passed as parameters are in degrees
        self.ra = ra
        self.dec = dec

        # ra and dec if passed as non Nne variables are converted into 
        # radians
        if self.dec is not None:
            self.dec = dec * np.pi / 180.0 
        if self.ra is not None:
            self.ra = ra * np.pi / 180.0

        self.lsstmwebv = EBVbase()
        self.ebvofMW = None
        if self.ra is not None and self.dec is not None:
            self.mwEBVfromMaps()
        return


    def set_MWebv(self, value):
        """
        if mwebv value is known, this can be used to set the attribute
        ebvofMW of the SNObject class to the value (float).

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

        .. note:: This function must be run after the class has attributes ra 
                  and dec set. In case it is run before this, the mwebv value 
                  will be set to None.

        """

        if self.ra is None or self.dec is None:
            return
        # else set skycoord
        self.skycoord = np.array([[self.ra], [self.dec]])
        self.ebvofMW  = self.lsstmwebv.calculateEbv(equatorialCoordinates=
                                                    self.skycoord)[0]
        return

    def SNObjectSED(self, time, wavelen=None, bandpassobject=None):
        '''
        return a sims_photutils.sed  object from the SN model at the requested
        time and wavelengths (or bandpassobject) with extinction
        from MW according to the SED extinction methods. If the sed is requested
        at times outside the validity range of the model, 0. is returned as the
        flux density. If the time is within the range of validity of the model,
        but the wavelength range requested is outside the range, then the
        returned fluxes are np.nan outside the range, and the model fluxes
        inside

        Parameters
        ----------
        time: float
            time of observation
        wavelen: `np.ndarray` of floats, optional, defaults to None
            array containing wavelengths in nm
        bandpassobject: `sims_photUtils.bandpass` object or list thereof,
                         optional, defaults to `None`

            if provided, overrides wavelen input and the SED is obtained at the
            wavelength values native to bandpass object.


        Returns
        -------
        `sims_photutils.sed` object containing the wavelengths and SED values
        from the SN at time time


        .. note: If both wavelen and bandpassobject are `None` then exception,
                 will be raised.


        Examples
        --------
        >>> sed = SN.SNObjectSED(time=0., wavelen=wavenm)
        >>> # Usual Sed methods and attributes can be accessed
        '''

        if wavelen is None and bandpassobject is None:
            raise ValueError('A non None input to either wavelen or bandpassobject must be provided')

        if bandpassobject is not None:
            if isinstance(bandpassobject, list):
                bp = bandpassobject[0]
            else:
                bp = bandpassobject
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

            wave = wavelen * 10.

            mask1 = wave > self.minwave()
            mask2 = wave < self.maxwave()
            mask = mask1 & mask2

            wave = wave[mask]


            flambda[mask] = self.flux(time=time, wave=wave) * 10.
            print flambda

        SEDfromSNcosmo = Sed(wavelen=wavelen, flambda=flambda)

        # Apply LSST extinction
        ax, bx = SEDfromSNcosmo.setupCCMab()
        if self.ebvofMW is None:
            raise ValueError('ebvofMW attribute cannot be None Type, and must\
                    be set using either ra, dec or by hand useing set_MWebv \
                    before this stage \n')

        SEDfromSNcosmo.addCCMDust(a_x=ax, b_x=bx, ebv=self.ebvofMW)
        return SEDfromSNcosmo
        
    def snObjectFlux(self, time, wave):
        '''
        returns the spectral flux density at the requested times and wavelengths

        '''
        newflux = np.empty(np.shape(wave), dtype=float)
        newflux[:] = np.nan
                    
        mask1 = wave > self.minwave()
        mask2 = wave < self.maxwave()
        mask = mask1 & mask2

        newflux[mask] = self.flux(time= time, wave=wave[mask])
        # newflux = np.where( newflux > 0., newflux, 0.) 
        return newflux

     
    def bandFluxes(self, bandpassobjects, time, phiarray=None):
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

        # CODE to call function returning the SED that is not working
        # Talk to Scott about it

        # SEDfromSNcosmo = self.SNObjectSED(bandpassobjects, time)

        # print '------------'
        # print type(SEDfromSNcosmo), SEDfromSNcosmo.__class__
        # print '------------'

        # self.parameters contains a list of values of SNModel as defined 
        # The parameters are z, t0, x1, c 
        # in sncosmo.Model. These can be set using the method `SNObject.set( ) 
        # inherited from sncosmo.Model

        # z = self.parameters[0] 
        #if self.parameters[0] > 1.2:
        #    return [np.nan]*len(bandpassobjects)

        filterwav = bandpassobjects[0].wavelen


        # print 'TIME ', time
        # print self.parameters
        
        # This is dE/dt/dA/dlambda 
        # Since lambda is expressed in nm
        sncosmofluxdensity = self.snObjectFlux(time=time, 
                                               wave=filterwav*10.)*10.
        SEDfromSNcosmo = Sed(wavelen=filterwav,
                             flambda=sncosmofluxdensity)

        # Apply LSST extinction
        ax, bx = SEDfromSNcosmo.setupCCMab()
        if self.ebvofMW is None:
            raise ValueError('ebvofMW attribute cannot be None Type, and must\
                    be set using either ra, dec or by hand useing set_MWebv \
                    before this stage \n')
        SEDfromSNcosmo.addCCMDust(a_x=ax, b_x=bx, ebv=self.ebvofMW)


        if phiarray is None:
            phiarray, dlambda = SEDfromSNcosmo.setupPhiArray(bandpassobjects)

        SEDfromSNcosmo.flambdaTofnu()
        
        return SEDfromSNcosmo.manyFluxCalc(phiarray, wavelen_step=wavelenstep)
