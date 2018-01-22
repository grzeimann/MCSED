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
import ssp
import cosmology
import emcee
import corner
import time

import numpy as np
import matplotlib.pyplot as plt


class Mcsed:
    def __init__(self, filter_matrix, ssp_spectra, ssp_ages, ssp_masses,
                 ssp_met, wave, sfh_name, dust_abs_name, data_fnu=None,
                 data_fnu_e=None, redshift=None, filter_flag=None,
                 input_spectrum=None, input_params=None, sigma_m=0.1,
                 nwalkers=40, nsteps=1000, true_fnu=None):
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
        wave : numpy array (1 dim)
            wavelength for SSP models and all model spectra
        sfh_class : class
            This is the input class for sfh.  Each class has a common attribute
            which is "sfh_class.nparams" for organizing the total model_params.
            Also, each class has a key function, sfh_class.evaluate(t), with
            the input of time in units of Gyrs
        dust_abs_class : class
            This is the input class for dust absorption.
        data_fnu : numpy array (1 dim)
            Photometry for data.  Length = (filter_flag == True).sum()
        data_fnu_e : numpy array (1 dim)
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
        self.ssp_met = ssp_met
        self.wave = wave
        self.sfh_class = getattr(sfh, sfh_name)()
        self.dust_abs_class = getattr(dust_abs, dust_abs_name)()
        self.ssp_class = getattr(ssp, 'fsps_freeparams')()
        self.SSP = None
        self.param_classes = ['sfh_class', 'dust_abs_class', 'ssp_class']
        self.data_fnu = data_fnu
        self.data_fnu_e = data_fnu_e
        self.redshift = redshift
        self.filter_flag = filter_flag
        self.input_spectrum = input_spectrum
        self.input_params = input_params
        self.sigma_m = sigma_m
        self.nwalkers = nwalkers
        self.nsteps = nsteps
        self.true_fnu = true_fnu
        if self.redshift is not None:
            self.set_new_redshift(self.redshift)

        # Set up logging
        self.setup_logging()

        # Time array for sfh
        self.age_eval = np.logspace(-3, 1, 4000)

    def set_new_redshift(self, redshift):
        ''' Setting redshift

        Parameters
        ----------
        redshift : float
            Redshift of the source for fitting
        '''
        self.redshift = redshift
        # Need luminosity distance to adjust ssp_spectra from 10pc to Dl
        self.Dl = cosmology.Cosmology().luminosity_distance(self.redshift)
        self.sfh_class.set_agelim(self.redshift)

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

    def get_filter_wavelengths(self):
        ''' FILL IN
        '''
        wave_avg = np.dot(self.wave, self.filter_matrix[:, self.filter_flag])
        return wave_avg

    def get_filter_fluxdensities(self):
        '''Convert a spectrum to photometric fluxes for a given filter set.
        The photometric fluxes will be in the same units as the spectrum.
        Ideally, the spectrum should be in microjanskies(lambda) such that
        the photometric fluxes will be in microjanskies.

        Returns
        -------
        f_nu : numpy array (1 dim)
            Photometric flux densities for an input spectrum
        '''
        f_nu = np.dot(self.spectrum, self.filter_matrix[:, self.filter_flag])
        return f_nu

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

        ######################################################################
        # SSP Parameters
        self.ssp_class.set_parameters_from_list(theta, start_value)
        start_value += self.ssp_class.nparams

    def get_ssp_spectrum(self):
        if self.ssp_class.fix_met:
            if self.SSP is not None:
                return self.SSP
        met_weights = np.zeros(self.ssp_met.shape)
        Z = np.log10(self.ssp_met)
        ind = np.searchsorted(Z, self.ssp_class.met)
        if ind == 0:
            met_weights[0] = 1.
        elif ind == len(Z):
            met_weights[ind-1] = 1.
        else:
            D = Z[ind] - Z[ind-1]
            met_weights[ind] = (self.ssp_class.met - Z[ind-1]) / D
            met_weights[ind-1] = 1. - met_weights[ind]
        self.SSP = np.dot(self.ssp_spectra, met_weights)
        return self.SSP

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
        # Collapse for metallicity
        SSP = self.get_ssp_spectrum()

        # Need star formation rate from observation back to formation
        sfr = self.sfh_class.evaluate(self.ssp_ages)
        ageval = 10**self.sfh_class.age

        # ageval sets limit on ssp_ages that are useable in model calculation
        sel = self.ssp_ages <= ageval

        # The weight is the time between ages of each SSP
        weight = np.diff(np.hstack([0, self.ssp_ages])) * 1e9 * sfr
        # Ages greater than ageval should have zero weight in csp
        weight[~sel] = 0

        # Cover the two cases where ssp_ages contains ageval and when not
        A = np.nonzero(self.ssp_ages <= ageval)[0][-1]
        select_too_old = np.nonzero(self.ssp_ages >= ageval)[0]
        if len(select_too_old):
            B = np.nonzero(self.ssp_ages >= ageval)[0][0]
            if A == B:
                spec_dustfree = np.dot(SSP, weight)
                mass = np.sum(weight * self.ssp_masses)
            else:
                lw = ageval - self.ssp_ages[A]
                wei = lw * 1e9 * np.interp(ageval, self.ssp_ages, sfr)
                weight[B] = wei
                spec_dustfree = np.dot(SSP, weight)
                mass = np.sum(weight * self.ssp_masses)
        else:
            spec_dustfree = np.dot(SSP, weight)
            mass = np.sum(weight * self.ssp_masses)

        # Need to correct for dust attenuation
        Alam = self.dust_abs_class.evaluate(self.wave)
        spec_dustobscured = spec_dustfree * 10**(-0.4 * Alam)

        # Redshift to observed frame
        csp = np.interp(self.wave, self.wave * (1. + self.redshift),
                        spec_dustobscured * (1. + self.redshift))
        # Correct spectra from 10pc to redshift of the source
        return csp / self.Dl**2, mass

    def lnprior(self):
        ''' Simple, uniform prior for input variables

        Returns
        -------
        0.0 if all parameters are in bounds, -np.inf if any are out of bounds
        '''
        flag = True
        for par_cl in self.param_classes:
            flag *= getattr(self, par_cl).prior()
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
        inv_sigma2 = 1.0 / (self.data_fnu_e**2 + (model_y * self.sigma_m)**2)
        chi2_term = -0.5 * np.sum((self.data_fnu - model_y)**2 * inv_sigma2)
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

    def get_init_walker_values(self, kind='ball', num=None):
        ''' Before running emcee, this function generates starting points
        for each walker in the MCMC process.

        Returns
        -------
        pos : np.array (2 dim)
            Two dimensional array with Nwalker x Ndim values
        '''
        # We need an initial guess for emcee so we take it from the model class
        # parameter values and deltas
        init_params, init_deltas, init_lims = [], [], []
        for par_cl in self.param_classes:
            init_params.append(getattr(self, par_cl).get_params())
            init_deltas.append(getattr(self, par_cl).get_param_deltas())
            if len(getattr(self, par_cl).get_param_lims()):
                init_lims.append(getattr(self, par_cl).get_param_lims())
        theta = list(np.hstack(init_params))
        thetae = list(np.hstack(init_deltas))
        theta_lims = np.vstack(init_lims)
        if num is None:
            num = self.nwalkers
        if kind == 'ball':
            pos = emcee.utils.sample_ball(theta, thetae, size=num)
        else:
            pos = (np.random.rand(num)[:, np.newaxis] *
                   (theta_lims[:, 1]-theta_lims[:, 0]) + theta_lims[:, 0])
        return pos

    def get_param_names(self):
        ''' Grab the names of the parameters for plotting

        Returns
        -------
        names : list
            list of all parameter names
        '''
        names = []
        for par_cl in self.param_classes:
            names.append(getattr(self, par_cl).get_names())
        names = list(np.hstack(names))
        return names

    def get_params(self):
        ''' Grab the the parameters in each class

        Returns
        -------
        vals : list
            list of all parameter values
        '''
        vals = []
        for par_cl in self.param_classes:
            vals.append(getattr(self, par_cl).get_params())
        vals = list(np.hstack(vals))
        return vals

    def get_param_lims(self):
        ''' Grab the limits of the parameters for making mock galaxies

        Returns
        -------
        limits : numpy array (2 dim)
            an array with parameters for rows and limits for columns
        '''
        limits = []
        for par_cl in self.param_classes:
            limits.append(getattr(self, par_cl).get_param_lims())
        limits = np.array(sum(limits, []))
        return limits

    def fit_model(self):
        ''' Using emcee to find parameter estimations for given set of
        data magnitudes and errors
        '''
        # Need to verify data parameters have been set since this is not
        # a necessity on initiation
        self.log.info('Fitting model using emcee')
        check_vars = ['data_fnu', 'data_fnu_e', 'redshift', 'filter_flag']
        for var in check_vars:
            if getattr(self, var) is None:
                self.error('The variable %s must be set first' % var)

        pos = self.get_init_walker_values(kind='ball')
        ndim = pos.shape[1]
        start = time.time()
        # Time to set up the sampler and run the mcmc
        sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.lnprob,
                                        a=2.0)

        # Do real run
        sampler.run_mcmc(pos, self.nsteps, rstate0=np.random.get_state())
        end = time.time()
        elapsed = end - start
        self.log.info("Total time taken: %0.2f s" % elapsed)
        self.log.info("Time taken per step per walker: %0.2f ms" %
                      (elapsed / (self.nsteps) * 1000. /
                       self.nwalkers))
        # Calculate how long the run should last
        tau = np.max(sampler.acor)
        burnin_step = int(tau*3)
        self.log.info("Mean acceptance fraction: %0.2f" %
                      (np.mean(sampler.acceptance_fraction)))
        self.log.info("AutoCorrelation Steps: %i, Number of Burn-in Steps: %i"
                      % (np.round(tau), burnin_step))
        new_chain = np.zeros((self.nwalkers, self.nsteps, ndim+2))
        new_chain[:, :, :-2] = sampler.chain
        self.chain = sampler.chain
        for i in xrange(len(sampler.blobs)):
            for j in xrange(len(sampler.blobs[0])):
                x = sampler.blobs[i][j]
                new_chain[j, i, -2] = np.where((np.isfinite(x)) * (x > 10.),
                                               np.log10(x), -99.)
        new_chain[:, :, -1] = sampler.lnprobability
        self.samples = new_chain[:, burnin_step:, :].reshape((-1, ndim+2))

    def spectrum_plot(self, ax, color=[0.996, 0.702, 0.031]):
        ''' Make spectum plot for current model '''
        spectrum, mass = self.build_csp()
        ax.plot(self.wave, spectrum, color=color, alpha=0.1)

    def add_subplots(self, fig, nsamples):
        ''' Add Subplots to Triangle plot below '''
        ax1 = fig.add_subplot(3, 1, 1)
        ax1.set_position([0.7, 0.80, 0.25, 0.15])
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_ylabel(r'SFR $M_{\odot} yr^{-1}$')
        ax1.set_xlabel('Lookback Time (Gyr)')
        ax1.set_xlim([10**-3,
                      10**self.sfh_class.age_lims[1]])
        ax1.set_ylim([1e-5, 1e3])
        ax2 = fig.add_subplot(3, 1, 2)
        ax2.set_position([0.7, 0.60, 0.25, 0.15])
        ax2.set_xscale('log')
        xtick_pos = [1000, 3000, 10000]
        xtick_lbl = ['1000', '3000', '10000']
        ax2.set_xticks(xtick_pos)
        ax2.set_xticklabels(xtick_lbl)
        ax2.set_xlim([1000, 20000])
        ax2.set_ylim([0, 8])
        ax2.set_ylabel(r'Dust Attenuation (mag)')
        ax2.set_xlabel(r'Wavelength $\AA$')
        ax3 = fig.add_subplot(3, 1, 3)
        ax3.set_position([0.38, 0.80, 0.25, 0.15])
        ax3.set_xscale('log')
        xtick_pos = [5000, 10000, 20000, 50000]
        xtick_lbl = ['0.5', '1', '2', '5']
        ax3.set_xticks(xtick_pos)
        ax3.set_xticklabels(xtick_lbl)
        ax3.set_xlim([3000, 80000])
        ax3.set_xlabel(r'Wavelength $\mu m$')
        ax3.set_ylabel(r'$F_{\nu}$ ($\mu$Jy)')
        rndsamples = 100
        for i in np.arange(rndsamples):
            ind = np.random.randint(0, nsamples.shape[0])
            self.set_class_parameters(nsamples[ind, :-2])
            self.sfh_class.plot(ax1)
            self.dust_abs_class.plot(ax2, self.wave)
            self.spectrum_plot(ax3)
        if self.input_params is not None:
            self.set_class_parameters(self.input_params)
            self.sfh_class.plot(ax1, color='k')
            self.dust_abs_class.plot(ax2, self.wave, color='k')
        wv = self.get_filter_wavelengths()
        if self.true_fnu is not None:
            p = ax3.scatter(wv, self.true_fnu, marker='o', s=25,
                            color=[0.216, 0.471, 0.749], zorder=9)
            p.set_facecolor('none')
        ax3.errorbar(wv, self.data_fnu, yerr=self.data_fnu_e, fmt='s',
                     fillstyle='none', markersize=2,
                     color=[0.510, 0.373, 0.529], zorder=10)
        ax3min = np.percentile(self.data_fnu, 5)
        ax3max = np.percentile(self.data_fnu, 95)
        ax3ran = ax3max - ax3min
        ax3.set_ylim([ax3min - 0.4 * ax3ran, ax3max + 0.4 * ax3ran])

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
        o = 0  # self.sfh_class.nparams
        names = self.get_param_names()[o:]
        names.append('Log Mass')
        if self.input_params is not None:
            truths = self.input_params[o:]
        else:
            truths = None
        percentilerange = [p for i, p in enumerate(self.get_param_lims())
                           if i >= o] + [[7, 11]]
        percentilerange = [.99] * len(names)
        fig = corner.corner(nsamples[:, o:-1], labels=names,
                            range=percentilerange,
                            truths=truths,
                            label_kwargs={"fontsize": 18}, show_titles=True,
                            title_kwargs={"fontsize": 16},
                            quantiles=[0.16, 0.5, 0.84], bins=30)
        # Adding subplots
        self.add_subplots(fig, nsamples)
        fig.set_size_inches(15.0, 15.0)
        fig.savefig("%s.png" % (outname), dpi=150)
        plt.close()

    def sample_plot(self, outname):
        ''' Make a sample plot

        Input
        -----
        outname : string
            The sample plot will be saved as "sample_{outname}.png"

        '''
        # Make selection for three sigma sample
        names = self.get_param_names()
        if self.input_params is not None:
            truths = self.input_params
        else:
            truths = None
        fig, ax = plt.subplots(self.chain.shape[2], 1, sharex=True,
                               figsize=(5, 2*self.chain.shape[2]))
        for i, a in enumerate(ax):
            for chain in self.chain[:, :, i]:
                a.plot(chain, 'k-', alpha=0.3)
            a.set_ylabel(names[i])
            if truths is not None:
                a.plot([0, self.chain.shape[1]], [truths[i], truths[i]], 'r--')
        fig.savefig("%s.png" % (outname))
        plt.close(fig)

    def add_fitinfo_to_table(self, percentiles, start_value=3, lnprobcut=7.5):
        chi2sel = (self.samples[:, -1] >
                   (np.max(self.samples[:, -1], axis=0) - lnprobcut))
        nsamples = self.samples[chi2sel, :-1]
        n = len(percentiles)
        for i, per in enumerate(percentiles):
            for j, v in enumerate(np.percentile(nsamples, per, axis=0)):
                self.table[-1][(i + start_value + j*n)] = v
