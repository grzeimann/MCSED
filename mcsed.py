""" SED fitting class using emcee for parameter estimation

    CURRENT LIMITATIONS:
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
                 input_spectrum=None, input_params=None, sigma_m=0.02):
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
        sigma_m : float
            Fractional error expected from the models.  This is used in
            the log likelihood calculation.  No model is perfect, and this is
            more or less a fixed parameter to encapsulate that.
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
        self.sigma_m = sigma_m

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
        The photometric fluxes will be in the same units as the spectrum.
        Ideally, the spectrum should be in microjanskies(lambda) such that
        the photometric fluxes will be in microjanskies.

        Returns
        -------
        mags : numpy array (1 dim)
            Photometric magnitudes for an input spectrum
        '''
        mags = np.dot(self.spectrum, self.filter_matrix[:, self.filter_flag])
        return mags

    def set_class_parameters(self, theta):
        ''' For a given set of model parameters, set the needed class variables
        related to SFH, dust attenuation, ect.

        Input
        -----
        theta : list
            list of input parameters for sfh, dust att., and dust em.
        '''
        start_value = 0
        ######################################################################
        # STAR FORMATION HISTORY
        self.sfh_class.set_parameters_from_list(theta, start_value)
        # Keeping track of theta index for age of model and other classes
        start_value += self.sfh_class.nparams

        ######################################################################
        # DUST ATTENUATION
        self.dust_abs_class.set_parameters_from_list(theta, start_value)
        start_value += self.dust_abs_class.nparams

    def build_csp(self):
        '''Build a composite stellar population model for a given star
        formation history, dust attenuation law, and dust emission law.

        Returns
        -------
        csp : numpy array (1 dim)
            Composite stellar population model
        mass : float
            Mass for csp given the SFH input
        '''
        sfh = self.sfh_class.evaluate(self.age_eval)
        ageval = self.sfh_class.age

        # ageval sets limit on ssp_ages that are useable in model calculation
        sel = self.ssp_ages <= ageval
        # Need star formation rate from observation back to formation
        sfr = np.interp(ageval - self.ssp_ages, self.age_eval, sfh)
        # The weight is the time between ages of each SSP
        weight = np.diff(np.hstack([0, self.ssp_ages])) * 1e9 * sfr
        # Ages greater than ageval should have zero weight in csp
        weight[~sel] = 0

        # Cover the two cases where ssp_ages contains ageval and when not
        A = np.nonzero(self.ssp_ages <= ageval)[0][-1]
        B = np.nonzero(self.ssp_ages >= ageval)[0][0]
        if A == B:
            spec_dustfree = np.dot(self.ssp_spectra, weight)
            mass = np.sum(weight * self.ssp_masses)
        else:
            lw = ageval - self.ssp_ages[A]
            wei = lw * 1e9 * np.interp(ageval, self.age_eval, sfh)
            weight[B] = wei
            spec_dustfree = np.dot(self.ssp_spectra, weight)
            mass = np.sum(weight * self.ssp_masses)

        # Need to correct for dust attenuation
        taulam = self.dust_abs_class.evaluate(self.wave)
        spec_dustobscured = spec_dustfree * np.exp(-1 * taulam)

        # Redshift to observed frame
        csp = np.interp(self.wave, self.wave * (1. + self.redshift),
                        spec_dustobscured)
        return csp, mass

    def lnprior(self):
        ''' Simple, uniform prior for input variables

        Returns
        -------
        0.0 if all parameters are in bounds, -np.inf if any are out of bounds
        '''
        sfh_flag = self.sfh_class.prior()
        dust_abs_flag = self.dust_abs_class.prior()
        flag = sfh_flag * dust_abs_flag
        if not flag:
            return -np.inf
        else:
            return 0.0

    def lnlike(self):
        ''' Calculate the log likelihood and return the value and stellar mass
        of the model

        Returns
        -------
        log likelihood, mass : float, float
            The log likelihood includes a chi2_term and a parameters term.
            The mass comes from building of the composite stellar population
        '''
        self.spectrum, mass = self.build_csp()
        model_y = self.get_filter_fluxdensities()
        inv_sigma2 = 1.0 / (self.data_magerrs**2 + (model_y * self.sigma_m)**2)
        chi2_term = -0.5 * np.sum((self.data_mags - model_y)**2 * inv_sigma2)
        parm_term = -0.5 * np.sum(np.log(1 / inv_sigma2))
        return (chi2_term + parm_term, mass)

    def lnprob(self, theta):
        ''' Calculate the log probabilty and return the value and stellar mass
        of the model

        Returns
        -------
        log prior + log likelihood, mass : float, float
            The log probability is just the sum of the logs of the prior and
            likelihood.  The mass comes from the building of the composite
            stellar population.
        '''
        self.set_class_parameters(theta)
        lp = self.lnprior()
        if np.isfinite(lp):
            lnl, mass = self.lnlike()
            return lp + lnl, mass
        else:
            return -np.inf, -np.inf
