""" MCSED - ssp.py

Single Stellar Population module for loading models

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""
import sfh
import sys

import numpy as np
import os.path as op
import scipy.interpolate as scint
from astropy.convolution import Gaussian1DKernel, convolve


def get_coarser_wavelength_fsps(wave, spec, redwave=1e5):
    '''
    smooth the spectrum with a gaussian kernel to improve 
    computational efficiency

    only affects the wavelength grid
    (the age and metallicity grids remain unchanged)

    Parameters
    ----------
    wave : numpy array (1d)
        initial wavelength grid
    spec : numpy array
        initial SSP grid over (wave, age, metallicity)
    redwave : float
        red wavelength cutoff (in Angstroms)
    '''
# WPBWPB: coarse_wavelength, how specified by user?
    sel = np.where((wave > 500) * (wave < redwave))[0]
    spec = spec[sel, :]
    wave = wave[sel]
    G = Gaussian1DKernel(25)
    nsel = np.where(np.abs(np.diff(wave)-0.9) < 0.5)[0]
    for i in np.arange(spec.shape[1]):
        spec[nsel, i] = convolve(spec[nsel, i], G)
    ndw = 12.
    nw = np.arange(wave[nsel[0]], wave[nsel[-1]+1], ndw)
    nwb = np.hstack([nw, wave[nsel[-1]+1]])
    nwave = np.hstack([wave[:nsel[0]], nw, wave[(nsel[-1]+1):]])
    nspec = np.zeros((len(nwave), spec.shape[1]))
    for i, sp in enumerate(spec.swapaxes(0, 1)):
        hsp = (np.histogram(wave[nsel], nwb, weights=sp[nsel])[0] /
               np.histogram(wave[nsel], nwb)[0])
        nspec[:, i] = np.hstack([sp[:nsel[0]], hsp, sp[(nsel[-1]+1):]])
    return nwave, nspec


def bin_ages_fsps(args, ages, spec):
    ''' FILL IN

    returns age (Gyr), blah
    '''
    sfh_class = getattr(sfh, args.sfh)()
    sel = ages >= 6.
    ages, spec = (ages[sel], spec[:, sel])
    weight = np.diff(np.hstack([0., 10**ages]))
    agebin = np.hstack([0., sfh_class.ages])
    nspec = np.zeros((spec.shape[0], len(sfh_class.ages)))
    for i in np.arange(len(sfh_class.ages)):
        sel = np.where((ages >= agebin[i]) * (ages < agebin[i+1]))[0]
        nspec[:, i] = np.dot(spec[:, sel], weight[sel]) / weight[sel].sum()
    return 10**(np.array(sfh_class.ages)-9.), nspec


def read_fsps_neb(filename):
    '''
    Returns
    -------
    Z : list (1 dim)
        metallicity grid in log solar units
    Age : list (1 dim)
        age grid in years (up to 10 Myr)
    logU : list (1 dim)
        ionization parameter grid
    spec : list (1 dim)
        elements are spectra (numpy array, 1 dim)
        relative line fluxes, normalized by number of ionizing photons
    wave : numpy array (1 dim)
        wavelength for each spectrum in Angstroms
    '''
    cnt = 0
    Z, Age, logU, spec = [], [], [], []
    with open(filename) as f:
        for lines in f:
            if cnt == 1:
                wave = np.array(lines.split(), dtype=float)
            if cnt > 1:
                l = lines.split()
                if len(l) == 3:
                    Z.append(float(l[0]))
                    Age.append(float(l[1]))
                    logU.append(float(l[2]))
                else:
                    spec.append(np.array(l, dtype=float))
            cnt += 1
    return Z, Age, logU, spec, wave


def read_fsps_ssp(filename):
    '''
    Returns
    -------
    ages : list (1 dim)
        ages in log years
    masses : list (1 dim)
    spec : list (2 dim)
        spectra in solar bolometric luminosity per Hz
    wave : numpy array (1 dim)
        wavelength for each spectrum spanning 90 Angstroms - 10000 microns
        in units of Angstroms
    '''
    cnt = 0
    ages, masses, spec = [], [], []
    with open(filename) as f:
        for lines in f:
            if cnt == 9:
                wave = np.array(lines.split(), dtype=float)
            if cnt > 9:
                l = lines.split()
                if len(l) == 4:
                    ages.append(float(l[0]))
                    masses.append(float(l[1]))
                else:
                    spec.append(np.array(l, dtype=float))
            cnt += 1
    return ages, masses, spec, wave


def read_fsps(args, metallicity):
    '''Read in the stellar population models from fsps for a given isochrone
    and metallicity.

    Parameters
    ----------
    args : class
        The args class from mcsed.parse_args()

    Returns
    -------
    ages : numpy array (1 dim)
        ages of each SPS (Gyr)
    masses : numpy array (1 dim)
        Remnant mass at a given age (solar masses)
    wave : numpy array (1 dim)
        wavelength grid for each spectrum in units of Angstroms
    spec : numpy array (2 dim)
        Spectra in f_nu (micro Janskies, i.e., 1e-29 ergs/s/cm^2/Hz) at 10pc
    '''
    pc10 = 10. * 3.08567758e18
    solar_microjansky = 3.826e33 * 1e29 / (4. * np.pi * pc10**2)
    filename = op.join('SSP', 'fsps_%s_%0.4f.spec' % (args.isochrone,
                                                      metallicity))
    if not op.exists(filename):
        print('Tried to open %s' % filename)
        print('Metallicity entered, %0.4f, does not match any of the %s '
              'isochrones of the %s models' % (args.metallicity,
                                               args.isochrone, args.ssp))
        print('Metallicity options [')
        for met in args.metallicity_dict[args.isochrone]:
            print('%0.4f ' % met)
        print(']')
        sys.exit(1)
    ages, masses, spec, wave = read_fsps_ssp(filename)
    spec = np.array(spec).swapaxes(0, 1) * solar_microjansky
    ages, masses = (np.array(ages), np.array(masses))
    # Total mass including remnants, so set to 1.
    sel = (ages <= 9.5) * (ages >= 6.)
    return 10**(ages[sel]-9), np.ones(ages[sel].shape), wave, spec[:, sel]

def get_nebular_emission(ages, wave, spec, logU, metallicity,
                         filename='nebular/ZAU_ND_pdva',
                         sollum=3.826e33, kind='both'):
    ''' 
    ages : numpy array (1 dim)
        ages of the SSP models
    wave : numpy array (1 dim)
        wavelength for SSP models
    spec : numpy array (3 dim)
        SSP spectrum for each age and each metallicity
    kind : str
        {'both', 'line', 'cont'}
        'line', 'cont' return only the line and continuum nebular emission
        'both' returns line and continuum nebular emission
    '''
    while kind not in ['line', 'cont', 'both']:
        kind = input("Invalid entry. Please enter 'line', 'cont', or 'both'")
    cont_file = filename + '.cont'
    lines_file = filename + '.lines'
    cont_res = [np.array(x) for x in read_fsps_neb(cont_file)]
    if kind != 'cont':
        lines_res = [np.array(x) for x in read_fsps_neb(lines_file)]
    # Make array of Z, age, U
    V = np.array([10**cont_res[0]*0.019, cont_res[1]/1e6,
                  cont_res[2]]).swapaxes(0, 1)
    if kind != 'line':
        # 1e48 factor avoids interpolating small numbers - removed later
        C = scint.LinearNDInterpolator(V, cont_res[3]*1e48)
    if kind != 'cont':
        L = scint.LinearNDInterpolator(V, lines_res[3]*1e48)
        garray = make_gaussian_emission(wave, lines_res[4])
# WPBWPB: uncomment
    nspec = spec * 0.
## WPBWPB delete: only when saving an array without gaussian emission:
##    print(nspec.shape)
#    nspec = nspec[0:len(lines_res[4]), :] 
#    garray = 1.
    for i, age in enumerate(ages):
# WPBWPB: if age of birth cloud is user-defined, adjust here:
        if age <= 1e-2:
            if kind != 'line':
                cont = C(metallicity, age*1e3, logU)
            if kind != 'cont':
                lines = L(metallicity, age*1e3, logU)
            qq = number_ionizing_photons(wave, spec[:, i]) / 1e48 * sollum
            if kind == 'both': 
                nspec[:, i] = (nspec[:, i] 
                               + np.interp(wave, cont_res[4], cont*qq)
                               + (garray * lines * qq).sum(axis=1))
            if kind == 'line':
## WPBWPB delete: only when saving an array without gaussian emission:
#                nspec[:, i] = lines * qq
# WPBWPB: uncomment - original code
                nspec[:, i] = (nspec[:, i] 
                               + (garray * lines * qq).sum(axis=1))
            if kind == 'cont':
                nspec[:, i] = (nspec[:, i] 
                               + np.interp(wave, cont_res[4], cont*qq))

## WPBWPB: only for case when no gaussian emission is included
#    np.savez('SSP_nongaussian', linewave0=lines_res[4], linespec0=nspec)
#    return 

    return nspec

def add_nebular_emission(ages, wave, spec, logU, metallicity,
                         filename='nebular/ZAU_ND_pdva',
                         sollum=3.826e33):
    cont_file = filename + '.cont'
    lines_file = filename + '.lines'
    cont_res = [np.array(x) for x in read_fsps_neb(cont_file)]
    lines_res = [np.array(x) for x in read_fsps_neb(lines_file)]
    # Make array of Z, age, U
    V = np.array([10**cont_res[0]*0.019, cont_res[1]/1e6,
                  cont_res[2]]).swapaxes(0, 1)
    C = scint.LinearNDInterpolator(V, cont_res[3]*1e48)
    L = scint.LinearNDInterpolator(V, lines_res[3]*1e48)
    garray = make_gaussian_emission(wave, lines_res[4])
    nspec = spec * 1.
    for i, age in enumerate(ages):
        if age <= 1e-2:
            cont = C(metallicity, age*1e3, logU)
            lines = L(metallicity, age*1e3, logU)
            qq = number_ionizing_photons(wave, spec[:, i]) / 1e48 * sollum
            nspec[:, i] = (nspec[:, i] + np.interp(wave, cont_res[4], cont*qq)
                           + (garray * lines * qq).sum(axis=1))
    return nspec

def collapse_emline_SSP(args, linewave, linespec, clight=2.99792e18):
    '''Speed up construction of emission line fluxes from the CSP

    Parameters
    ----------
    args : dictionary
        user-passed arguments from config.py and command line
    linewave : numpy array (1 dim)
        wavelength for SSP models
    linespec : numpy array (3 dim)
        SSP spectrum for each age and each metallicity

    Returns
    -------
    linewave : numpy array (1 dim)
        collapsed wavelength grid of emission lines
    linespec : numpy array (3 dim)
        collapsed SSP spectrum for each age and each metallicity
        in units of ergs / s / cm^2 at 10 pc
    ''' 

    if not args.use_emline_flux:
        return np.zeros(2), linespec[0:2,:,:]

# loop through all emission line spectra for all ages, metallicities
# WPBWPB: generalize such that does not assume only grid over ages and metallicities, but maybe ionization parameter (or arbitrary number of properties)
# maybe I want to raise an error?
    emlines = args.emline_list_dict.keys() 
    emwaves = np.array(args.emline_list_dict.values())[:,0]
## WPBWPB delete
#    print('here are emwaves, type: %s, %s' % (emwaves, type(emwaves)))
    ssp_emline_collapsed = linespec[0:len(emlines),:,:]
## WPBWPB delete
#    waves_collapsed = []
    # loop through emission line spectra for all ages, metallicities
    dims = linespec.shape
    for i in range(dims[1]): # age array
        # if no emission lines at this age, skip metallicity grid
        if np.max(linespec[:,i,:]) <= 0:
            empty = np.zeros( (len(emlines), dims[2]) )
            ssp_emline_collapsed[:,i,:] = empty
            continue
        for j in range(dims[2]): # metallicity array
            spec = linespec[:,i,j]
            lineflux = []
            for emline in emlines:
## WPBWPB delete
#                print('here"s the emline dict and your element:')
#                print(args.emline_list_dict)
#                print([type(val) for val in args.emline_list_dict.values()])
#                print(emline)
#                print(args.emline_list_dict[emline])
#                print(type(args.emline_list_dict[emline]))

                w = args.emline_list_dict[emline][0]
                indx = np.searchsorted(linewave, w)
                if spec[indx] <= 0:
                    lineflux.append(0.)
                    continue
                # find indices when target line goes to zero
                indx_zero = np.where(spec==0)[0]
                indx_lo = indx - np.min( indx - indx_zero[indx_zero<indx])
                indx_hi = indx + np.min( indx_zero[indx_zero>indx] - indx)
                dnu = np.diff( clight / linewave[indx_lo:indx_hi+2])
                y = spec[indx_lo:indx_hi+1]
                # convert flux density in micro-Jy at 10 pc 
                # to total flux (ergs/s/cm2) at 10pc
                lineflux.append(np.dot( y, np.abs(dnu) ) / 1e29)
            ssp_emline_collapsed[:,i,j] = lineflux

    # return modified emission line wavelength and SSP grid: 
    # fluxes at discrete wavelengths
    return np.array(emwaves), ssp_emline_collapsed


def make_gaussian_emission(wavebig, wave, stddev=1., clight=2.99792e18):
    gspec = np.zeros((len(wavebig), len(wave)))
    G = Gaussian1DKernel(stddev).array
    mid = len(G) / 2
    dw = np.diff(wavebig)
    for i, w in enumerate(wave):
        xl = np.argmin(np.abs(wavebig - w))
        if (xl > mid) and ((xl+mid) < len(wavebig)):
            gspec[(xl-mid):(xl+mid+1), i] = G / clight * w**2 / dw[xl]
    return gspec


def number_ionizing_photons(wave, spectrum, clight=2.99792e18,
                            hplanck=6.626e-27):
    nu = clight / wave
    xlim = np.searchsorted(wave, 912., side='right')
    x1 = np.abs(np.diff(nu[:xlim])) # dnu
    x2 = spectrum[:xlim] / nu[:xlim] # micro-Jy at 10 pc / nu
    x3 = (x2[:-1] + x2[1:]) / 2.
    # WPBWPB add to comments
    # the sum of x1 * x3 gives micro-Jy at 10 pc / nu
    # dividing by hplanck gives (number of ionizing photons * 1e-29 / cm2 at 10 pc)
    return np.sum(x1 * x3) / hplanck


def read_ssp(args):
    ''' Read in SPS model and return ages, masses, wavelength, and spectra

    Parameters
    ----------
    args : class
        The args class from mcsed.parse_args()
WPBWPB: need to determine units of all outputs
WPBWPB: operate under assumption that spec, linespec are in same units
    '''
    import matplotlib.pyplot as plt
    plt.ioff()

    s, ls, m = ([], [], [])
    for met in args.metallicity_dict[args.isochrone]:
        if args.ssp.lower() == 'fsps':
            ages, masses, wave, spec = read_fsps(args, met)
        # carry emission lines only, for measuring line fluxes
# WPBWPB: only carry linespec if going to measure emission lines?
        if args.add_nebular:
            linespec = get_nebular_emission(ages, wave, spec, args.logU,
                                        met, kind='line')
            spec = add_nebular_emission(ages, wave, spec, args.logU,
                                        met)
        else:
            linespec = np.zeros(np.shape(spec))

# WPBWPB add comment
        if args.sfh == 'empirical' or args.sfh == 'empirical_direct':
            ages0 = ages.copy()
            ages, spec = bin_ages_fsps(args, np.log10(ages0)+9., spec)
            ages9, linespec = bin_ages_fsps(args, np.log10(ages0)+9., linespec)
        masses = np.ones(ages.shape)
        # do not smooth the emission line grid
        wave0 = wave.copy()
        if args.fit_dust_em:
            wave, spec = get_coarser_wavelength_fsps(wave0, spec, redwave=200e4)
        else:
            wave, spec = get_coarser_wavelength_fsps(wave0, spec)
        wei = (np.diff(np.hstack([0., ages])) *
               getattr(sfh, args.sfh)().evaluate(ages))
        s.append(spec)
        ls.append(linespec)
# WPBWPB following line appears unused... 
        m.append(np.dot(spec, wei)/spec.shape[1])

    # for plotting purposes only
    fig = plt.figure(figsize=(8, 8))
    import seaborn as sns
    colors = sns.color_palette("coolwarm", s[-5].shape[1])
    wei = np.diff(np.hstack([0., ages]))
    for i in np.arange(s[-5].shape[1]):
        plt.plot(wave, s[-5][:, i] * wei[i] / 1e8, color=colors[i])
    plt.xlim([900., 40000.])
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1e-5, 10])
    plt.savefig('template_spectra_plot.png')
    plt.close(fig)

    spec = np.moveaxis(np.array(s), 0, 2)
    linespec = np.moveaxis(np.array(ls), 0, 2)
    metallicities = args.metallicity_dict[args.isochrone]

## WPBWPB delete
#    linespec0 = linespec.copy()

    # Collapse the emission line SSP grid
    linewave, linespec = collapse_emline_SSP(args, wave0, linespec) 

## WPBWPB delete
#    np.savez('SSP', wave=wave, spec=spec, linewave=linewave, linespec=linespec, linewave0=wave0, linespec0=linespec0)

## WPB delete
#    print('these are ages: %s' % ages)

    return ages, masses, wave, spec, np.array(metallicities), linewave,linespec


class fsps_freeparams:
    ''' Allowing metallicity to be free -0.39'''
# WPBWPB why is there an -0.39 in here...? might appear elsewhere too - grep it
    def __init__(self, met=-0.7, met_lims=[-1.98, 0.2], met_delta=0.3,
                 fix_met=False):
        ''' Initialize this class

        Parameters
        ----------
        TODO Fill these in
        '''
        self.met = met
        self.met_lims = met_lims
        self.met_delta = met_delta
        self.fix_met = fix_met

    def get_nparams(self):
        ''' Return number of parameters '''
        if self.fix_met:
            return 0
        else:
            return 1

    def get_params(self):
        ''' Return current parameters '''
        l = []
        if not self.fix_met:
            l.append(self.met)
        return l

    def get_param_lims(self):
        ''' Return current parameters limits '''
        l = []
        if not self.fix_met:
            l.append(self.met_lims)
        return l

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        l = []
        if not self.fix_met:
            l.append(self.met_delta)
        return l

    def get_names(self):
        ''' Return names of each parameter '''
        l = []
        if not self.fix_met:
            l.append('Log Z')
        return l

    def prior(self):
        ''' Uniform prior based on boundaries '''
        flag = (self.met >= self.met_lims[0]) * (self.met <= self.met_lims[1])
        return flag

    def set_parameters_from_list(self, input_list, start_value):
        ''' Set parameters from a list and a start_value

        Parameters
        ----------
        input_list : list
            list of input parameters (could be much larger than number of
            parameters to be set)
        start_value : int
            initial index from list to read out parameters
        '''
        if not self.fix_met:
            self.met = input_list[start_value]
