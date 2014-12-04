"""
Class describing the SN object itself. The SN object derives from SNCosmo.Model and has additional properties such 
as ra, dec.  
It also has additional methods :
-     calc_mags which use the magnitude calculations in LSST stack
-     extinction which use the extinction calculations in LSST stack

- Usage:  (after setups described in the readme.rst) python snObject.py 

"""
import sncosmo
import numpy as np
import os
import matplotlib.pyplot as plt

from astropy.units import Unit
from astropy.coordinates import SkyCoord
from sncosmo import Model

from lsst.sims.photUtils.Photometry import PhotometryStars, Sed, Bandpass

dustmaproot = os.getenv('SIMS_DUSTMAPS_DIR')
map_dir = os.path.join(dustmaproot, 'DustMaps')

wavelenstep = 0.1
plot=False


def getlsstbandpassobjs(loadsncosmo=True, loadcatsim=True):
    bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
    banddir = os.path.join(os.getenv('THROUGHPUTS_DIR'), 'baseline')
    lsstbands = []
    lsstbp = {}

    #wavelenstep = 0.1
    for band in bandPassList:
        # setup sncosmo bandpasses
        bandfname = banddir + "/total_" + band  + '.dat'
        if loadsncosmo:
            # Usually the next two lines can be merged, but there is an astropy bug currently.
            numpyband = np.loadtxt(bandfname)
            sncosmoband = sncosmo.Bandpass(wave=numpyband[:,0], trans=numpyband[:,1], wave_unit=Unit('nm'), 
                                           name='LSST'+band)
            sncosmo.registry.register(sncosmoband, force=True)
        if loadcatsim:                     
            # Now load LSST bandpasses for catsim
            lsstbp[band] = Bandpass()
            lsstbp[band].readThroughput(bandfname, wavelen_step=wavelenstep)
            lsstbands.append(lsstbp[band])
    fifilterfigs, filterax = plt.subplots()
    if plot:
        for band in bandPassList:
            b = sncosmo.get_bandpass('LSST' + band)
            filterax.plot (b.wave, b.trans, '-k', lw=2.0)
            filterax.set_xlabel(r'$\lambda$, ($\AA$)')
            filterax.set_ylabel(r'transmission')
        plt.show()

    if loadcatsim:
        return lsstbands 
    else:
        return None


class SNObject (Model):

    # def __init__(self, ra , dec):
    def __init__(self, ra=None, dec=None):
        Model.__init__(self, source="salt2-extended", effects=[sncosmo.CCM89Dust()], effect_names=['mw'],
                effect_frames=['obs'])
        # default value of mwebv = 0. ie. no extinction 
        self.set(mwebv=0.)
        self._ra = ra
        self._dec = dec
        if ra is None or dec is None:
            pass
        else:
            self.setmwebv()
        #skycoords = SkyCoord(ra, dec, unit='deg')
        #self._mwebv = sncosmo.get_ebv_from_map(skycoords, mapdir=map_dir, interpolate=False)
        #self.snIaModel.set(**params)
        self._seed = None
        #    self._z = None
        return


    @property
    def seed(self):
        return self._seed


    @property
    def ra(self):
        #self._ra = ra
        return self._ra 

    @property
    def dec(self):
        #dec = self._dec
        return self._dec

    def setmwebv(self):
        ra = self._ra
        dec = self._dec
        skycoords = SkyCoord(ra, dec, unit='deg')
        self._mwebv = sncosmo.get_ebv_from_map(skycoords, mapdir=map_dir, interpolate=False)
        return

    def get_SNparams(self):

        hundredyear = 100*365.0
        np.random.seed(self.seed)
        t0 = np.random.uniform(-hundredyear/2.0, hundredyear/2.0)
        c = np.random.normal(0., 0.3)
        x1 = np.random.normal(0., 1.0)
        mabs = np.random.normal(-19.3, 0.3)
        self.snIaModel.set(z=self.z, c=c, x1=x1)
        self.snIaModel.set_source_peakabsmag(mabs, 'bessellb', 'ab')
        x0 = self.snIaModel.get('x0')

        return np.array([t0, x0, x1, c])

    def lsstbandmags(self, lsstbands, time):

        filterwav = lsstbands[0].wavelen
        print "within lsstbandmags", len(filterwav)
        #print 'lsstbandmags: ', time, filterwav
        SEDfromSNcosmo = Sed(wavelen=filterwav, flambda=self.flux(time=time, wave=filterwav*10.)*10.)
        SEDfromSNcosmo.synchronizeSED(wavelen_min=filterwav[0], wavelen_max=filterwav[-2], wavelen_step=wavelenstep)
        print "SED :", len(SEDfromSNcosmo.wavelen)
        ax, bx = SEDfromSNcosmo.setupCCMab()
        print len(ax), len(bx), len(SEDfromSNcosmo.flambda), len(SEDfromSNcosmo.fnu)
        print "********************", self._mwebv
        SEDfromSNcosmo.addCCMDust(a_x=ax, b_x=bx, ebv=self._mwebv)
        phiarray, dlambda = SEDfromSNcosmo.setupPhiArray(lsstbands)
        SEDfromSNcosmo.synchronizeSED(wavelen_min=filterwav[0], wavelen_max=filterwav[-2], wavelen_step=wavelenstep)

        return SEDfromSNcosmo.manyMagCalc(phiarray, wavelen_step=wavelenstep) 



if __name__ == "__main__":
    """ 
    Example code for writing a light curve in multiple bands

    """

    import numpy as np
    import matplotlib.pyplot as plt

    ra = 204.
    dec = -30.
    SN = SNObject(ra, dec)
    SN.set(x0=1.847e-6, x1=0.1, c=0., z =0.2)
    #print SN
    SNCosmoSN = SNObject(ra, dec)
    SNCosmoSN.set(x0=1.847e-6, x1=0.1, z=0.2, mwebv=SN._mwebv)
    lsstbands = getlsstbandpassobjs()
    sncosmobands = ['LSSTu', 'LSSTg', 'LSSTr', 'LSSTi', 'LSSTz']
    w  = sncosmo.get_bandpass(sncosmobands[0]).wave
    l = [] 
    for time in np.arange(-20.,50.,1.0):
        if time < 24. and time > 17.:
            plt.plot(w, SN.flux(time=time, wave=w))
            plt.yscale('log')
            #plt.show()
        t = time*np.ones(len(sncosmobands))
        t.tolist()
        #print time
        #print SN.flux(time=time, wave = [3000., 5000., 12000.])#, SN.flux(time=time, wave=12000.)
        #try:
        #    x = SN.lsstbandmags(lsstbands, time=time)
        #except:
        #    print "Something wierd"
        #    print time
        #    SN.flux(wave=sncosmobands.get_bandpass().wave, time=time, zpsys='ab')
        #print "passed ?"
        y = SNCosmoSN.bandmag(band=['LSSTg','LSSTr'], time=[time,time], magsys='ab')
        #print "yes !"
        e =  [time ]
        e += x.tolist() 
        e += y.tolist()
        l.append(e)

    np.savetxt('lc.dat', np.array(l))
