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
import sfh
import dust_abs
import cosmology
import emcee
import triangle
import time

import numpy as np
import matplotlib.pyplot as plt


class Mcsed:
    def __init__(self, filter_matrix, ssp_spectra, ssp_ages, ssp_masses,
                 wavelength, sfh_name, dust_abs_name, data_mags=None,
                 data_magerrs=None, redshift=None, filter_flag=None,
                 input_spectrum=None, input_params=None, sigma_m=0.02,
                 nwalkers=40, nsteps=500):
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
        nwalkers : int
            The number of walkers for emcee when fitting a model
        nsteps : int
            The number of steps each walker will make when fitting a model
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
        self.nwalkers = nwalkers
        self.nsteps = nsteps
        self.Dl = cosmology.Cosmology().luminosity_distance(self.redshift)

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
            Composite stellar population model at self.redshift
        mass : float
            Mass for csp given the SFH input
        '''
        # TODO include D(z)**2 modification to ssp_spectra
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
        # Correct spectra from 10pc to redshift of the source
        return csp / (4. * np.pi * self.Dl**2), mass

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

    def get_init_walker_values(self):
        ''' Before running emcee, this function generates starting points
        for each walker in the MCMC process.

        Returns
        -------
        pos : np.array (2 dim)
            Two dimensional array with Nwalker x Ndim values
        '''
        # We need an initial guess for emcee so we take it from the model class
        # parameter values and deltas
        init_params = []
        init_deltas = []
        param_classes = ['sfh_class', 'dust_abs_class']
        for par_cl in param_classes:
            init_params.append(getattr(self, par_cl).get_params())
            init_deltas.append(getattr(self, par_cl).get_param_deltas())
        theta = list(np.hstack(init_params))
        thetae = list(np.hstack(init_deltas))
        pos = emcee.utils.sample_ball(theta, thetae, size=self.nwalkers)
        return pos

    def get_param_names(self):
        ''' Grab the names of the parameters for plotting

        Returns
        -------
        names : list
            list of all parameter names
        '''
        names = []
        param_classes = ['sfh_class', 'dust_abs_class']
        for par_cl in param_classes:
            names.append(getattr(self, par_cl).get_names())
        names = list(np.hstack(names))
        return names

    def fit_model(self):
        ''' Using emcee to find parameter estimations for given set of
        data magnitudes and errors
        '''
        # Need to verify data parameters have been set since this is not
        # a necessity on initiation
        self.log.info('Fitting model using emcee')
        check_vars = ['data_mags', 'data_magerrs', 'redshift', 'filter_flag']
        for var in check_vars:
            if getattr(self, var) is None:
                self.error('The variable %s must be set first' % var)

        pos = self.get_init_walker_values()
        ndim = pos.shape[1]
        start = time.time()
        # Time to set up the sampler and run the mcmc
        sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.lnprob,
                                        a=1.2)
        initial_steps = 100
        sampler.run_mcmc(pos, initial_steps)
        maxprobsearch = sampler.lnprobability == np.max(sampler.lnprobability)
        indx = np.nonzero(maxprobsearch)[0][0]
        indy = np.nonzero(maxprobsearch)[1][0]
        theta = sampler.chain[indx, indy, :]
        self.set_class_parameters(theta)
        pos = self.get_init_walker_values()

        # Reset the chain to remove the burn-in samples.
        sampler.reset()

        # Do real run
        sampler.run_mcmc(pos, self.nsteps, rstate0=np.random.get_state())
        end = time.time()
        elapsed = end - start
        self.log.info("Total time taken: %0.2f s" % elapsed)
        self.log.info("Time taken per step: %0.2f ms" %
                      (elapsed / (initial_steps + self.nsteps)))
        # Calculate how long the run should last
        tau = np.max(sampler.acor)
        burnin_step = np.round(tau*5)
        self.log.info("Mean acceptance fraction: %0.2f" %
                      (np.mean(sampler.acceptance_fraction)))
        self.log.info("AutoCorrelation Steps: %i, Number of Burn-in Steps: %i"
                      % (np.round(tau), burnin_step))
        new_chain = np.zeros((self.nwalkers, self.nsteps, ndim+2))
        new_chain[:, :, :-2] = sampler.chain
        for i in xrange(len(sampler.blobs)):
            for j in xrange(len(sampler.blobs[0])):
                new_chain[j, i, -2] = np.log10(sampler.blobs[i][j])
        new_chain[:, :, -1] = sampler.lnprobability
        self.samples = new_chain[:, burnin_step:, :].reshape((-1, ndim+2))

    def triangle_plot(self, outname, lnprobcut=7.5):
        ''' Make a triangle corner plot for samples from fit

        Input
        -----
        outname : string
            The triangle plot will be saved as "triangle_{outname}.png"
        lnprobcut : float
            Some of the emcee chains include outliers.  This value serves as
            a cut in log probability space with respect to the maximum
            probability.  For reference, a Gaussian 1-sigma is 2.5 in log prob
            space.
        '''
        # Make selection for three sigma sample
        chi2sel = (self.samples[:, -1] >
                   (np.max(self.samples[:, -1], axis=0) - lnprobcut))
        nsamples = self.samples[chi2sel, :]
        names = self.get_param_names()
        names.append('Log Mass')
        percentilerange = [0.99] * len(names)
        fig = triangle.corner(nsamples[:, :-1], labels=names,
                              range=percentilerange, truths=self.input_params,
                              label_kwargs={"fontsize": 18}, show_titles=True,
                              title_kwargs={"fontsize": 16},
                              quantiles=[0.16, 0.5, 0.84], bins=25)
        fig.set_size_inches(15.0, 15.0)
        fig.savefig("triangle_%s.png" % (outname), dpi=150)
        plt.close()
