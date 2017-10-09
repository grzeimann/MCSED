""" MCSED


1) CURRENT LIMITATIONS:
       A) Constant metallicity for input SSP
       B) Dust Emission is ad hoc from Draine and Li (2007)
   OPTIONAL FITTED PARAMETERS:
       A) SFH
           a) tau_sfh, age, a, b, c
       B) Dust law
           b) tau_dust, delta, Eb
   OUTPUT PRODUCTS:
       A) XXX Plot

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import logging
import numpy as np
import sfh
import dust_abs

__all__ = ["Mcsed"]


class Mcsed:
    def __init__(self, filter_matrix, ssp_spectra, ssp_ages, ssp_masses,
                 wavelength, sfh_name, dust_abs_name, data_mags=None,
                 data_magerrs=None, redshift=None, filter_flag=None,
                 input_spectrum=None, input_params=None):
        ''' Initialize the Mcsed class.

        Init
        ----
        filter_matrix : numpy array (2 dim)
            The filter_matrix has rows of wavelength and columns for each
            filter (can be much larger than the filters used for fitting)
        ssp_spectra : numpy array (2 dim)
            single stellar population spectrum for each age in ssp_ages
        ssp_ages : numpy array (1 dim)
            ages of the SSP models
        ssp_masses : numpy array (1 dim)
            remnant masses of the SSP models
        wavelength : numpy array (1 dim)
            wavelength for SSP models and all model spectra
        sfh_class : class
            This is the input class for sfh.  Each class has a common attribute
            which is "sfh_class.nparams" for organizing the total model_params.
            Also, each class has a key function, sfh_class.evaluate(t), with
            the input of time in units of Gyrs
        dust_abs_class : class
            This is the input class for dust absorption.
        data_mags : numpy array (1 dim)
            Photometry for data.  Length = (filter_flag == True).sum()
        data_magerrs : numpy array (1 dim)
            Photometric errors for data
        redshift : float
            Redshift of the source
        filter_flag : numpy array (1 dim)
            Length = filter_matrix.shape[1], True for filters matching data
        input_spectrum : numpy array (1 dim)
            F_nu(wave) for input
        input_params : list
            input parameters for modeling.  Intended for testing fitting
            procedure.
        '''
        # Initialize all argument inputs
        self.filter_matrix = filter_matrix
        self.ssp_spectra = ssp_spectra
        self.ssp_ages = ssp_ages
        self.ssp_masses = ssp_masses
        self.wavelength = wavelength
        self.sfh_class = getattr(sfh, sfh_name)()
        self.dust_abs_class = getattr(dust_abs, dust_abs_name)()
        self.data_mags = data_mags
        self.data_magerrs = data_magerrs
        self.redshift = redshift
        self.filter_flag = filter_flag
        self.input_spectrum = input_spectrum
        self.input_params = input_params

        # Set up logging
        self.setup_logging()

        # Time array for sfh
        self.age_eval = np.logspace(-3, 2, 2000)

    def setup_logging(self):
        '''Setup Logging for MCSED

        Builds
        -------
        self.log : class
            self.log.info() is for general print and self.log.error() is
            for raise cases
        '''
        self.log = logging.getLogger('mcsed')
        if not len(self.log.handlers):
            # Set format for logger
            fmt = '[%(levelname)s - %(asctime)s] %(message)s'
            fmt = logging.Formatter(fmt)
            # Set level of logging
            level = logging.INFO
            # Set handler for logging
            handler = logging.StreamHandler()
            handler.setFormatter(fmt)
            handler.setLevel(level)
            # Build log with name, mcsed
            self.log = logging.getLogger('mcsed')
            self.log.setLevel(logging.DEBUG)
            self.log.addHandler(handler)

    def get_filter_fluxdensities(self):
        '''Convert a spectrum to photometric fluxes for a given filter set.
        The photometric fluxes with be in the same units as the spectrum.
        Ideally, the spectrum should be in microjanskies(lambda) such that
        the photometric fluxes with be in microjanskies.

        Returns
        -------
        mags : numpy array (1 dim)
            Photometric magnitudes for an input spectrum
        '''
        mags = np.dot(self.spectrum, self.filter_matrix[:, self.filter_flag])
        return mags

    def build_csp(self, theta):
        '''Build a composite stellar population model for a given star
        formation history, dust attenuation law, and dust emission law.

        Input
        -----
        theta : list
            list of input parameters for sfh, dust att., and dust em.

        Returns
        -------
        csp : numpy array (1 dim)
            Composite stellar population model
        mass : float
            Mass for csp given the SFH input
        '''
        start_value = 0
        # Set sfh_class parameters from theta and keep track of theta index
        self.sfh_class.set_parameters_from_list(theta, start_value)
        start_value += self.sfh_class.nparams
        # Evaluate sfh_class
        sfh = self.sfh_class.evaluate(self.age_eval)
        # Get age of model and keep track of theta index
        ageval = 10**(theta[start_value])
        start_value += 1
        # Take only ages < age of galaxy for modeling
        sel = self.ages <= ageval
        # SFR interpolated from the finer sfh
        sfr = np.interp(ageval - self.ages, self.age_eval, sfh)
        # The weight is time between ages of each SSP
        weight = np.diff(np.hstack([0, self.ages])) * 1e9 * sfr
        # ages > ageval have zero weight
        weight[~sel] = 0
        # Cover the two cases where there is a fractional age range or not
        A = np.nonzero(self.ages <= ageval)[0][-1]
        B = np.nonzero(self.ages >= ageval)[0][0]
        if A == B:
            spec_dustfree = np.dot(self.ssp_spectra, weight)
            mass = np.sum(weight * self.masses)
        else:
            lw = ageval - self.ages[A]
            wei = lw * 1e9 * np.interp(ageval, self.age_eval, sfh)
            weight[B] = wei
            spec_dustfree = np.dot(self.ssp_spectra, weight)
            mass = np.sum(weight * self.masses)
        # Set dust attenuation parameters from theta
        self.dust_abs_class.set_parameters_from_list(theta, start_value)
        start_value += self.dust_abs_class.nparams
        # Evaluate dust_abs_class
        taulam = self.dust_abs_class.evaluate(self.wave)
        spec_dustobscured = spec_dustfree * np.exp(-1 * taulam)
        # Redshift
        csp = np.interp(self.wave, self.wave * (1. + self.redshift),
                        spec_dustobscured)
        return csp, mass

    def testbounds(theta):
        return False
    
    
    def lnprior(theta):
        flag = testbounds(theta)
        if flag:
            return -np.inf
        else:
            return 0.0
    
    
    def lnlike(theta, ages, seds, masses, wave, zobs, y, yerr, filters, sigma_m):
        flag = testbounds(theta)
        if flag:
            return -np.inf, -np.inf
        else:
            spec, mass = build_csp(theta, ages, seds, masses, wave, zobs)
            model_y = get_filter_fluxdensities(spec, filter_flag, filter_matrix)
            inv_sigma2 = 1.0/(yerr**2+(model_y*sigma_m)**2)
            return (-0.5*np.sum((y-model_y)**2*inv_sigma2) 
                        - 0.5*np.sum(np.log(1/inv_sigma2)), mass)
    
    
    def lnprob(theta, ages, seds, masses, wave, zobs, y, yerr, filters):
        lp = lnprior(theta)
        lnl , mass = lnlike(theta, ages, seds, masses, wave, y, yerr, zobs, 
                            filters)
        if not np.isfinite(lp):
            return -np.inf, -np.inf
        return lp+lnl, mass