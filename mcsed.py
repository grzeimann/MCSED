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
import dust_emission
import ssp
import cosmology
import emcee
import corner
import time
# WPBWPB delete astrpy table
from astropy.table import Table

import numpy as np
import matplotlib.pyplot as plt
plt.ioff()

#WPBWPB re organize the arguments (aesthetic purposes)
class Mcsed:
    def __init__(self, filter_matrix, ssp_spectra,
                 emlinewave, ssp_emline, ssp_ages, ssp_masses,
                 ssp_met, wave, sfh_class, dust_abs_class, dust_em_class,
                 data_fnu=None, data_fnu_e=None, 
                 data_emline=None, data_emline_e=None, emline_dict=None,
                 redshift=None,
                 filter_flag=None, input_spectrum=None, input_params=None,
                 sigma_m=0.1, nwalkers=40, nsteps=1000, true_fnu=None):
        ''' Initialize the Mcsed class.

        Init
        ----
        filter_matrix : numpy array (2 dim)
            The filter_matrix has rows of wavelength and columns for each
            filter (can be much larger than the filters used for fitting)
        ssp_spectra : numpy array (3 dim)
            single stellar population spectrum for each age in ssp_ages
            and each metallicity in ssp_met 
        emlinewave : numpy array (1 dim)
            Rest-frame wavelengths of requested emission lines (emline_dict)
            Corresponds to ssp_emline
        ssp_emline : numpy array (3 dim)
            Emission line SSP grid spanning emlinewave, age, metallicity
            Only includes requested emission lines (from emline_dict)
            Only used for calculating model emission line strengths
            Spectral units are ergs / s / cm2 at 10 pc
        ssp_ages : numpy array (1 dim)
            ages of the SSP models
        ssp_masses : numpy array (1 dim)
            remnant masses of the SSP models
        ssp_met : numpy array (1 dim)
            metallicities of the SSP models
        wave : numpy array (1 dim)
            wavelength for SSP models and all model spectra
        sfh_class : str
            Converted from str to class in initialization
            This is the input class for sfh.  Each class has a common attribute
            which is "sfh_class.nparams" for organizing the total model_params.
            Also, each class has a key function, sfh_class.evaluate(t), with
            the input of time in units of Gyrs
        dust_abs_class : str 
            Converted from str to class in initialization
            This is the input class for dust absorption.
        dust_em_class : str
            Converted from str to class in initialization
            This is the input class for dust absorption.
        data_fnu : numpy array (1 dim)
            Photometry for data.  Length = (filter_flag == True).sum()
WPBWPB units + are dimensions correct??
        data_fnu_e : numpy array (1 dim)
            Photometric errors for data
        data_emline : Astropy Table (1 dim)
            Emission line fluxes in units ergs / cm2 / s
        data_emline_e : Astropy Table (1 dim)
            Emission line errors in units ergs / cm2 / s
        emline_dict : dictionary
            Keys are emission line names (str)
            Values are rest-frame wavelength in Angstroms (float)
        use_emline_flux : bool
            If emline_dict contains emission lines, set to True. Else, False
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
        true_fnu : WPBWPB FILL IN
        '''
        # Initialize all argument inputs
        self.filter_matrix = filter_matrix
        self.ssp_spectra = ssp_spectra
        self.emlinewave = emlinewave
        self.ssp_emline = ssp_emline
        self.ssp_ages = ssp_ages
        self.ssp_masses = ssp_masses
        self.ssp_met = ssp_met
        self.wave = wave
        self.dnu = np.abs(np.hstack([0., np.diff(2.99792e18 / self.wave)]))
        self.sfh_class = getattr(sfh, sfh_class)()
        self.dust_abs_class = getattr(dust_abs, dust_abs_class)()
        self.ssp_class = getattr(ssp, 'fsps_freeparams')()
        self.dust_em_class = getattr(dust_emission, dust_em_class)()
# WPBWPB: describe SSP, lineSSP in comments... 
# ssp_spectra span many metallicities, SSP only span ages
        self.SSP = None
        self.lineSSP = None
# WPBWPB is ssp_class still a thing?
        self.param_classes = ['sfh_class', 'dust_abs_class', 'ssp_class',
                              'dust_em_class']
        self.data_fnu = data_fnu
        self.data_fnu_e = data_fnu_e
        self.data_emline = data_emline
        self.data_emline_e = data_emline_e
        self.emline_dict = emline_dict
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

    def remove_waverange_filters(self, wave1, wave2, restframe=True):
        '''Remove filters in a given wavelength range

        Parameters
        ----------
        wave1 : float
            start wavelength of masked range (in Angstroms)
        wave2 : float
            end wavelength of masked range (in Angstroms)
        restframe : bool
            if True, wave1 and wave2 correspond to rest-frame wavelengths
        '''
        wave1, wave2 = np.sort([wave1, wave2])
        if restframe:
            wave_factor = 1. + self.redshift
        else:
            wave_factor = 1.
        loc1 = np.searchsorted(self.wave, wave1 * wave_factor)
        loc2 = np.searchsorted(self.wave, wave2 * wave_factor)
        # account for the case where indices are the same
        if (loc1 == loc2):
            loc2+=1
        maxima = np.max(self.filter_matrix, axis=0)
        try:
            newflag = np.max(self.filter_matrix[loc1:loc2, :], axis=0) < maxima * 0.1
        except ValueError:
            return
        maximas = np.max(self.filter_matrix[:, self.filter_flag], axis=0)
        newflags = np.max(self.filter_matrix[loc1:loc2, self.filter_flag], axis=0) < maximas * 0.1
        self.filter_flag = self.filter_flag * newflag
        if self.true_fnu is not None:
            self.true_fnu = self.true_fnu[newflags]
        self.data_fnu = self.data_fnu[newflags]
        self.data_fnu_e = self.data_fnu_e[newflags]


    def get_filter_wavelengths(self):
        ''' FILL IN
        '''
        wave_avg = np.dot(self.wave, self.filter_matrix[:, self.filter_flag])
        return wave_avg

    def get_filter_fluxdensities(self):
        '''Convert a spectrum to photometric fluxes for a given filter set.
        The photometric fluxes will be in the same units as the spectrum.
        The spectrum is in microjanskies(lambda) such that
        the photometric fluxes will be in microjanskies.

        Returns
        -------
        f_nu : numpy array (1 dim)
            Photometric flux densities for an input spectrum
        '''
# WPBWPB delete
        print('shape of spectrum, filter_matrix, filter_flag:')
        print((self.spectrum.shape, self.filter_matrix.shape, self.filter_flag.shape))
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
# WPBWPB modify: pass a dust_abs_birthcloud keyword, see if its the same, blah

        ######################################################################
        # SSP Parameters
        self.ssp_class.set_parameters_from_list(theta, start_value)
        start_value += self.ssp_class.nparams

        ######################################################################
        # DUST EMISSION
        self.dust_em_class.set_parameters_from_list(theta, start_value)
        start_value += self.dust_em_class.nparams

    def get_ssp_spectrum(self):
        '''
        Calculate SSP for an arbitrary metallicity (self.ssp_class.met) given a
        model grid for a range of metallicities (self.ssp_met)

        Returns
        -------
        SSP : 2-d array
            Single stellar population models for each age in self.ages
        '''
        if self.ssp_class.fix_met:
            if self.SSP is not None:
## WPBWPB delete
#                print('self.SSP is not None!')
                return self.SSP, self.lineSSP
        Z = np.log10(self.ssp_met)
        z = self.ssp_class.met + np.log10(0.019)
        X = Z - z
        wei = np.exp(-(X)**2 / (2. * 0.15**2))
        wei /= wei.sum()
        self.SSP = np.dot(self.ssp_spectra, wei)
# WPBWPB: this is where I would relax logU criteria, same as metallicity
# careful: needs to be dealt with when measuring line fluxes originally

        # only treat the emission line grid if requested
## WPBWPB delete
#        print('shape of emline SSP before/after wei')
#        print(self.ssp_emline.shape)
        if self.use_emline_flux:
            self.lineSSP = np.dot(self.ssp_emline, wei)
        else:
            self.lineSSP = self.ssp_emline[:,:,0]
## WPBWPB delete
#        print(self.lineSSP.shape)
        return self.SSP, self.lineSSP

    def build_dustfree_CSP(self, sfr, ageval, age_birth):
        '''WPBWPB FILL IN

        Parameters
        ----------
        sfr : WPBWPB units and form 
        ageval : float
            Current age of the object in Gyr
        age_birth : float
            Longevity of a birth cloud in Gyr

        Returns
        -------
        spec_dustfree : 
        spec_birth_dustfree :
        linespec_dustfree : 
        mass :        
        '''
        # ageval sets limit on ssp_ages that are useable in model calculation
        # age_birth separates birth cloud and diffuse components
# WPBWPB delete -- ageval, ssp_ages, age_birth are in units Gyr
        sel = (self.ssp_ages > age_birth) & (self.ssp_ages <= ageval)
        sel_birth = (self.ssp_ages <= age_birth) & (self.ssp_ages <= ageval)
        sel_age = self.ssp_ages <= ageval

        # The weight is the time between ages of each SSP
        weight = np.diff(np.hstack([0, self.ssp_ages])) * 1e9 * sfr
## WPBWPB delete
#        weight_orig = weight.copy()
        weight_birth = weight.copy()
        weight_age = weight.copy()
        # Ages greater than ageval should have zero weight in CSP
        # weight should only include populations younger than ageval
        # and older than age_birth
        # weight_birth should only include populations younger than ageval
        # and no older than age_birth
        # weight_age only considers the age of the system (for mass)
        weight[~sel] = 0
        weight_birth[~sel_birth] = 0
        weight_age[~sel_age] = 0

        # Cover the two cases where ssp_ages contains ageval and when not
        # A: index of last acceptable SSP age
        A = np.nonzero(self.ssp_ages <= ageval)[0][-1]
        # indices of SSP ages that are too old
        select_too_old = np.nonzero(self.ssp_ages >= ageval)[0]
        if len(select_too_old):
            # B: index of first SSP that is too old
            B = select_too_old[0]
            # only adjust weight if ageval falls between two SSP age gridpoints
            if A != B:
                lw = ageval - self.ssp_ages[A]
                wei = lw * 1e9 * np.interp(ageval, self.ssp_ages, sfr)
                if ageval > age_birth:
                    weight[B] = wei
                if ageval <= age_birth:
                    weight_birth[B] = wei
                weight_age[B] = wei

        # Adjust weights of the young component
        # Cover two cases where ssp_ages contains age_birth and when not
        # A: index of last acceptable SSP age
        A = np.nonzero(self.ssp_ages <= age_birth)[0][-1]
        # indices of SSP ages that are too old
        select_too_old = np.nonzero(self.ssp_ages >= age_birth)[0]
        if (len(select_too_old)>0): # & (ageval>=age_birth):
            # B: index of first SSP that is too old
            B = select_too_old[0]
            if A != B:
                lw = age_birth - self.ssp_ages[A]
                wei = lw * 1e9 * np.interp(age_birth, self.ssp_ages, sfr)
                if ageval > age_birth:
                    weight[B] = weight_age[B] - wei
                if ageval >= age_birth:
                    weight_birth[B] = wei
                else:
                    weight_birth[B] = weight_age[B]

### WPBWPB delete
##        print('this is ageval, age birth:   %s,  %s' % (ageval, age_birth))
#        # A summary table...
#        t=Table()
#        t['ageval'] = [ageval]*len(weight)
#        t['age_birth'] = [age_birth]*len(weight)
#        t['ssp_ages'] = self.ssp_ages
#        t['weight_orig'] = weight_orig
#        t['weight_young'] = weight_birth
#        t['weight_old'] = weight
#        t['weight_age'] = weight_age
#        t.write('CSP_weights_ageval%s_birth%s.dat' % (ageval, age_birth),format='ascii') 
#        return

        # Finally, do the matrix multiplication
        spec_dustfree = np.dot(self.SSP, weight)
        spec_birth_dustfree = np.dot(self.SSP, weight_birth)
        linespec_dustfree = np.dot(self.lineSSP, weight_birth)
        mass = np.sum(weight_age * self.ssp_masses)

        return spec_dustfree, spec_birth_dustfree, linespec_dustfree, mass


    def build_csp(self, sfr=None):
        '''Build a composite stellar population model for a given star
        formation history, dust attenuation law, and dust emission law.

        In addition to the returns it also modifies a lineflux dictionary

        Returns
        -------
        csp : numpy array (1 dim)
            Composite stellar population model at self.redshift
WPBWPB units??
        mass : float
            Mass for csp given the SFH input
        '''
        # Collapse for metallicity
        SSP, lineSSP = self.get_ssp_spectrum()

        # Need star formation rate from observation back to formation
        if sfr is None:
            sfr = self.sfh_class.evaluate(self.ssp_ages)
        ageval = 10**self.sfh_class.age # Gyr

        # Treat the birth cloud and diffuse component separately
# WPBWPB: may want to modify: have this as user-defined setting...
        age_birth = 10**-2 # Gyr 

### WPBWPB delete - both in Gyr
##        print('this is the ageval: %s' % ageval)
##        print('this is ssp ages: %s' % self.ssp_ages)
##        return

## WPBWPB delete
##        for ageval, age_birth in [ [0.008, 0.011 ], [0.005011872336272725, 0.011 ], [0.01, 0.011 ], [0.14, 0.011 ], [0.0116, 0.011 ], [0.1, 0.011 ], [0.0145, 0.01 ], [0.011, 0.011 ], [0.01, 0.01 ] ]:
#        for ageval, age_birth in [ [0.01116, 0.011] ]:
#            self.build_dustfree_CSP(sfr, ageval, age_birth)
#        return

        # Get dust-free CSPs, properly accounting for ages
        dustfree_CSP = self.build_dustfree_CSP(sfr, ageval, age_birth)
        spec_dustfree, spec_birth_dustfree, linespec_dustfree, mass = dustfree_CSP 

        # Need to correct spectrum for dust attenuation
        Alam = self.dust_abs_class.evaluate(self.wave)
        spec_dustobscured = spec_dustfree * 10**(-0.4 * Alam)

        # Correct the corresponding birth cloud spectrum separately
# WPBWPB: check which law using for the birth cloud
# if attenuating it directly tied to overall dust law,
# modified by coefficient between EBV_stars ~ gas, get it here
        Alam_birth = Alam / self.dust_abs_class.EBV_stars_gas
        spec_birth_dustobscured = spec_birth_dustfree * 10**(-0.4 * Alam_birth)

        # Combine the young and old components
        spec_dustfree += spec_birth_dustfree
        spec_dustobscured += spec_birth_dustobscured

        # recompute attenuation for observed wavelength of emission lines
        emwaves = self.emlinewave * (1. + self.redshift)
        Alam_emline = (self.dust_abs_class.evaluate(emwaves,new_wave=True)
                       / self.dust_abs_class.EBV_stars_gas)
## WPBWPB delete
#        print(emwaves)
#        print(Alam_emline)
#        print('shape of linespec_dustfree, emwaves : (%s, %s)' % (linespec_dustfree.shape, emwaves.shape))
        linespec_dustobscured = linespec_dustfree * 10**(-0.4*Alam_emline)
# else, use a screen model

# WPB: exclude dust emission component altogether? Does it make a difference?
        # Change in bolometric Luminosity
        L_bol = (np.dot(self.dnu, spec_dustfree) -
                 np.dot(self.dnu, spec_dustobscured))

        # Add dust emission
        spec_dustobscured += L_bol * self.dust_em_class.evaluate(self.wave)

        # Redshift to observed frame
        csp = np.interp(self.wave, self.wave * (1. + self.redshift),
                        spec_dustobscured * (1. + self.redshift))

        # Update dictionary of modeled emission line fluxes (observed)
        linefluxCSPdict = {}
        if self.use_emline_flux:
            for emline in self.emline_dict.keys():
                indx = np.argmin(np.abs(self.emlinewave 
                                        - self.emline_dict[emline]))
                flux = linespec_dustobscured[indx]
                # Correct flux from 10pc to redshift of source
                pc10cm = 10. * 3.08567758e18
                Dlcm = self.Dl * pc10cm
                linefluxCSPdict[emline] = linespec_dustobscured[indx] / Dlcm**2
        self.linefluxCSPdict = linefluxCSPdict

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

        # compare input and model emission line fluxes
        emline_term = 0.0
        emline_weight = 1.
        if self.use_emline_flux:
            # if all lines have null line strengths, ignore 
            if not min(self.data_emline) == max(self.data_emline) == -99:
## WPBWPB delete
#                print('this is emline_dict:' + str(self.emline_dict))
#                print('this is emline data, error:')
#                print(self.data_emline)
#                print(self.data_emline_e)
                for emline in self.emline_dict.keys():
                    if self.data_emline['%s_FLUX' % emline] > -99: # null value
                        emline_wave = self.emline_dict[emline]
                        model_lineflux = self.linefluxCSPdict[emline] 
                        lineflux  = self.data_emline['%s_FLUX' % emline]
                        elineflux = self.data_emline_e['%s_ERR' % emline]
                        emline_term += (-0.5 * (model_lineflux - lineflux)**2 /
                                        elineflux**2.) * emline_weight
## WPBWPB: delete
#                print('this is emline and term:')
#                print(emline)
#                print(emline_term)

        model_y = self.get_filter_fluxdensities()
        inv_sigma2 = 1.0 / (self.data_fnu_e**2 + (model_y * self.sigma_m)**2)
        chi2_term = -0.5 * np.sum((self.data_fnu - model_y)**2 * inv_sigma2)
        parm_term = -0.5 * np.sum(np.log(1 / inv_sigma2))
        return (chi2_term + parm_term + emline_term, mass)

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

    def get_derived_params(self):
        ''' These are not free parameters in the model, but are instead
        calculated from free parameters
        '''

        t20 = None
        t50 = None
        sfr10 = None
        sfr100 = None


    def spectrum_plot(self, ax, color=[0.996, 0.702, 0.031], alpha=0.1):
        ''' Make spectum plot for current model '''
        self.spectrum, mass = self.build_csp()
        ax.plot(self.wave, self.spectrum, color=color, alpha=alpha)

    def add_sfr_plot(self, ax1):
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_ylabel(r'SFR $M_{\odot} yr^{-1}$')
        ax1.set_xlabel('Lookback Time (Gyr)')
        ax1.set_xticks([1e-3, 1e-2, 1e-1, 1])
        ax1.set_xticklabels(['1 Myr', '10 Myr', '100 Myr', '1 Gyr'])
        ax1.set_yticks([1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3])
        ax1.set_yticklabels(['0.001', '0.01', '0.1', '1', '10', '100', '1000'])
        ax1.set_xlim([10**-3, 10**self.sfh_class.age_lims[1]])
        ax1.set_ylim([1e-3, 1e3])

    def add_dust_plot(self, ax2):
        ax2.set_xscale('log')
        xtick_pos = [1000, 3000, 10000]
        xtick_lbl = ['1000', '3000', '10000']
        ax2.set_xticks(xtick_pos)
        ax2.set_xticklabels(xtick_lbl)
        ax2.set_xlim([1000, 20000])
        ax2.set_ylim([0, 8])
        ax2.set_ylabel(r'Dust Attenuation (mag)')
        ax2.set_xlabel(r'Wavelength $\AA$')

    def add_spec_plot(self, ax3):
# WPBWPB: adjust wavelength range, depending on whether dust emission is fit
        ax3.set_xscale('log')
        if self.dust_em_class.fixed:
            xtick_pos = [3000, 5000, 10000, 20000, 40000]
            xtick_lbl = ['0.3', '0.5', '1', '2', '4']
            xlims = [3000, 50000]
        else:
            xtick_pos = [3000, 5000, 10000, 40000, 100000, 1000000]
            xtick_lbl = ['0.3', '0.5', '1', '4', '10', '100']
            xlims = [3000, 2000000]
        ax3.set_xticks(xtick_pos)
        ax3.set_xticklabels(xtick_lbl)
        ax3.set_xlim(xlims)
        ax3.set_xlabel(r'Wavelength $\mu m$')
        ax3.set_ylabel(r'$F_{\nu}$ ($\mu$Jy)')

    def add_subplots(self, ax1, ax2, ax3, nsamples):
        ''' Add Subplots to Triangle plot below '''
        wv = self.get_filter_wavelengths()
        rndsamples = 200
# WPBWPB edit: might not need hbm list anymore
        sp, fn, hbm = ([], [], [])
        for i in np.arange(rndsamples):
            ind = np.random.randint(0, nsamples.shape[0])
            self.set_class_parameters(nsamples[ind, :-2])
            self.sfh_class.plot(ax1, alpha=0.1)
            self.dust_abs_class.plot(ax2, self.wave, alpha=0.1)
            self.spectrum_plot(ax3, alpha=0.1)
            fnu = self.get_filter_fluxdensities()
            sp.append(self.spectrum * 1.)
            fn.append(fnu * 1.)
# WPB edit: plotting HBeta line
# used to have self.hbflux = self.measure_hb() --> changed
#            hbm.append(self.hbflux * 1.)
        # Plotting median value:
        self.medianspec = np.median(np.array(sp), axis=0)
#        self.hbmedian = np.median(hbm)
        ax3.plot(self.wave, self.medianspec, color='dimgray')
        self.fluxwv = wv
        self.fluxfn = np.median(np.array(fn), axis=0)
        ax3.scatter(self.fluxwv, self.fluxfn, marker='x', s=200,
                    color='dimgray', zorder=8)
        chi2 = (1. / (len(self.data_fnu) - 1) *
                (((self.data_fnu - self.fluxfn) / self.data_fnu_e)**2).sum())
        if self.input_params is not None:
            self.set_class_parameters(self.input_params)
            self.sfh_class.plot(ax1, color='k', alpha=1.0)
            self.dust_abs_class.plot(ax2, self.wave, color='k', alpha=1.0)
            self.spectrum_plot(ax3, color='k', alpha=0.5)
        if self.true_fnu is not None:
            p = ax3.scatter(wv, self.true_fnu, marker='o', s=150,
                            color=[0.216, 0.471, 0.749], zorder=9)
            p.set_facecolor('none')
        ax3.errorbar(wv, self.data_fnu, yerr=self.data_fnu_e, fmt='s',
                     fillstyle='none', markersize=15,
                     color=[0.510, 0.373, 0.529], zorder=10)
        
        sel = np.where((wv > 3000.) * (wv < 50000.))[0]
        ax3min = np.percentile(self.data_fnu[sel], 5)
        ax3max = np.percentile(self.data_fnu[sel], 95)
        ax3ran = ax3max - ax3min
        ax3.set_ylim([ax3min - 0.4 * ax3ran, ax3max + 0.4 * ax3ran])
        ax3.text(4200, ax3max + 0.2 * ax3ran, r'${\chi}_{\nu}^2 = $%0.2f' % chi2)

    def triangle_plot(self, outname, lnprobcut=7.5, imgtype='png'):
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
        imgtype : string
            The file extension of the output plot
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
        percentilerange = [.95] * len(names)
        fig = corner.corner(nsamples[:, o:-1], labels=names,
                            range=percentilerange,
                            truths=truths, truth_color='gainsboro',
                            label_kwargs={"fontsize": 18}, show_titles=True,
                            title_kwargs={"fontsize": 16},
                            quantiles=[0.16, 0.5, 0.84], bins=30)
        # Adding subplots
        ax1 = fig.add_subplot(3, 1, 1)
        ax1.set_position([0.7, 0.60, 0.25, 0.15])
        ax2 = fig.add_subplot(3, 1, 2)
        ax2.set_position([0.7, 0.40, 0.25, 0.15])
        ax3 = fig.add_subplot(3, 1, 3)
        ax3.set_position([0.38, 0.80, 0.57, 0.15])
        self.add_sfr_plot(ax1)
        self.add_dust_plot(ax2)
        self.add_spec_plot(ax3)
        self.add_subplots(ax1, ax2, ax3, nsamples)
# WPB edit: printing HBeta line flux on the figure
# used to have self.hbflux = self.measure_hb() --> changed
#        if self.sfh_class.hblim is not None:
#            fig.text(.5, .75, r'H$\beta$ input: %0.2f' %
#                     (self.sfh_class.hblim * 1e17), fontsize=18)
#        fig.text(.5, .70, r'H$\beta$ model: %0.2f' % (self.hbmedian * 1e17),
#                 fontsize=18)
        # fig.set_size_inches(15.0, 15.0)
        fig.savefig("%s.%s" % (outname, imgtype), dpi=150)
        plt.close(fig)

    def sample_plot(self, outname, imgtype='png'):
        ''' Make a sample plot

        Input
        -----
        outname : string
            The sample plot will be saved as "sample_{outname}.png"
        imgtype : string
            The file extension of the output plot

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
            if i == len(ax)-1:
                a.set_xlabel("Step")
        fig.savefig("%s.%s" % (outname, imgtype))
        plt.close(fig)

    def add_fitinfo_to_table(self, percentiles, start_value=3, lnprobcut=7.5):
        ''' Assumes that "Ln Prob" is the last column in self.samples
        '''
        chi2sel = (self.samples[:, -1] >
                   (np.max(self.samples[:, -1], axis=0) - lnprobcut))
        nsamples = self.samples[chi2sel, :-1]
        n = len(percentiles)
        for i, per in enumerate(percentiles):
            for j, v in enumerate(np.percentile(nsamples, per, axis=0)):
                self.table[-1][(i + start_value + j*n)] = v
        return (i + start_value + j*n)

    def add_truth_to_table(self, truth, start_value):
        for i, tr in enumerate(truth):
            self.table[-1][start_value + i + 1] = tr
