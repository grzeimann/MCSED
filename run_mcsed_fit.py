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

import argparse as ap
import numpy as np
import os.path as op
import logging
import config
from dust_abs import noll
# from dust_em import draine_li
from sfh import double_powerlaw
from ssp import read_ssp
from astropy.io import fits
from mcsed import Mcsed


def setup_logging():
    '''Setup Logging for MCSED, which allows us to track status of calls and
    when errors/warnings occur.

    Returns
    -------
    log : class
        log.info() is for general print and log.error() is for raise cases
    '''
    log = logging.getLogger('mcsed')
    if not len(log.handlers):
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
        log = logging.getLogger('mcsed')
        # FIXME find out if this is needed
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)
    return log


def parse_args(argv=None):
    '''Parse arguments from commandline or a manually passed list

    Parameters
    ----------
    argv : list
        list of strings such as ['-f', 'input_file.txt', '-s', 'default.ssp']

    Returns
    -------
    args : class
        args class has attributes of each input, i.e., args.filename
        as well as astributes from the config file
    '''

    parser = ap.ArgumentParser(description="MCSED",
                               formatter_class=ap.RawTextHelpFormatter)

    parser.add_argument("-f", "--filename",
                        help='''File to be read for galaxy data''',
                        type=str, default=None)

    parser.add_argument("-s", "--ssp",
                        help='''SSP Models, default fsps''',
                        type=str, default=None)

    parser.add_argument("-z", "--metallicity",
                        help='''Metallicity for SSP models, 0.02 is solar''',
                        type=float, default=None)

    parser.add_argument("-i", "--isochrone",
                        help='''Isochrone for SSP model, e.g. padova''',
                        type=str, default=None)

    parser.add_argument("-t", "--test",
                        help='''Test script with fake data''',
                        action="count", default=0)

    args = parser.parse_args(args=argv)

    # Use config values if none are set in the input
    arg_inputs = ['ssp', 'metallicity']
    for arg_i in arg_inputs:
        if getattr(args, arg_i) is None:
            setattr(args, arg_i, getattr(config, arg_i))

    # Copy list of config values to the args class
    config_copy_list = ['metallicity_dict', 'mock_masses', 'mock_redshift',
                        'mock_dust_tau', 'mock_dust_delta', 'mock_dust_bump',
                        'mock_sfh_tau', 'mock_sfh_b', 'mock_sfh_c',
                        'filt_dict', 'catalog_filt_dict',
                        'filter_matrix_name']

    for con_copy in config_copy_list:
        setattr(args, con_copy, getattr(config, con_copy))

    args.log = setup_logging()

    return args


def build_filter_matrix(args, wave):
    '''Build a filter matrix with each row being an index of wave and
    each column being a unique filter.  This makes computation from spectra
    to magnitudes quick and easy.

    Parameters
    ----------
    args : class
        The args class is carried from function to function with information
        from command line input and config.py
    wave : numpy array
        The wave array corresponds to the wavelengths of the SSP models being
        used.

    Returns
    -------
    Fil_matrix : numpy array (2 dim)
        As mentioned above, the Fil_matrix has rows of wavelength and
        columns for each filter in args.filt_dict/config.filt_dict
    '''
    if op.exists(args.filter_matrix_name):
        return np.loadtxt(args.filter_matrix_name)
    else:
        nfilters = len(args.filt_dict)
        Fil_matrix = np.zeros((len(wave), nfilters))
        for i in np.arange(nfilters):
            wv, through = np.loadtxt('FILTERS', args.filt_dict[i], unpack=True)
            new_through = np.interp(wave, wv, through, 0.0, 0.0)
            Fil_matrix[:, i] = new_through/np.sum(new_through)
        np.savetxt(args.filter_matrix_name, Fil_matrix)
        return Fil_matrix


def get_test_filters(args):
    '''Used in test mode, this function loops through args.filt_dict and sets
    a flag to true if the filter is in args.test_filter_dict or false if it
    is not.  This filter_flag is used later in the quick calculation of
    filter magnitudes.

    Parameters
    ----------
    args : class
        The args class is carried from function to function with information
        from command line input and config.py

    Returns
    -------
    filter_flag : numpy array (bool)
        Explained above.
    '''
    nfilters = len(args.filt_dict)
    filter_flag = np.zeros((nfilters,), dtype=bool)
    for i in args.filt_dict.keys():
        if i in args.test_filter_dict:
            filter_flag[i] = True
    return filter_flag


def read_input_file(args):
    '''This function reads a very specific input file and joins it with
    archived 3dhst catalogs.  The input file should have the following columns:
    FIELD, ID, Z

    Parameters
    ----------
    args : class
        The args class is carried from function to function with information
        from command line input and config.py

    Returns
    -------
    y : numpy array (2 dim)
        Photometric magnitudes from the 3DHST survey for each input source
    yerr : numpy array (2 dim)
        Photometric errors in magnitudes
    z : numpy array (1 dim)
        Redshift from the file returned as a numpy array
    flag : numpy array (2 dim)
        Flag set to True for filters in the catalog_filter_dict in config.py
    '''
    F = np.loadtxt(args.filename)
    nobj = F.shape[0]
    fields = ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']
    name_base = '_3dhst.v4.1.cat.FITS'
    field_dict = {}
    for field in fields:
        field_dict[field] = fits.open(op.join('3dhst_catalogs',
                                              field+name_base))[1]
    nfilters = len(args.filt_dict)
    y = np.zeros((nobj, nfilters))
    yerr = np.zeros((nobj, nfilters))
    flag = np.zeros((nobj, nfilters), dtype=bool)
    z = F[:, 2]
    # convert from mag_zp = 25 to microjanskies (mag_zp = 23.9)
    fac = 10**(-0.4*(25.0-23.9))
    for i, datum in enumerate(F):
        loc = datum[0].lower()
        for j, ind in enumerate(args.filt_dict.keys()):
            colname = "f_"+args.catalog_filter_dict[ind]
            ecolname = "e_"+args.catalog_filter_dict[ind]
            if colname in field_dict[loc].columns.names:
                y[i, j] = field_dict[loc].data[colname][int(datum[1])]*fac
                yerr[i, j] = field_dict[loc].data[ecolname][int(datum[1])]*fac
                flag[i, j] = True
            else:
                y[i, j] = 0.0
                yerr[i, j] = 0.0
                flag[i, j] = False
    return y, yerr, z, flag


def get_filter_fluxdensities(spec, filter_flag, filter_matrix):
    '''Convert a spectrum to photometric fluxes for a given filter set

    Parameters
    ----------
    spec : numpy array (2 dim)
        The spectrum of a model to be converted in magnitudes
    filter_flag : numpy array (1 dim)
        List of boolean the length of the number of filter in args.filt_dict
        True values are for filters in the catalog.
    filter_matrix : numpy array (2 dim)
        The filter_matrix has rows of wavelength and
        columns for each filter in args.filt_dict/config.filt_dict

    Returns
    -------
    mags : numpy array (1 dim)
        Photometric magnitudes for an input spectrum
    '''
    mags = np.dot(spec, filter_matrix[:, filter_flag])
    return mags


def build_csp(theta, ages, seds, masses, wave, zobs, sfh_func=double_powerlaw):
    '''Build a composite stellar population model for a set of input parameters

    Parameters
    ----------
    theta : list
        SED parameter list for SFH and dust absorption
    ages : numpy array (1 dim)
        Time in Gyr for the SSP models
    seds : numpy array (2 dim)
        SSP models for different ages and rows of wavelength
    masses : numpy array (1 dim)
        Remnant mass of each SSP models for each age
    wave : numpy array (1 dim)
        Wavelength for all SSP model
    zobs : float
        Redshift for model source

    Returns
    -------
    csp : numpy array (1 dim)
        Composite stellar population model
    mass : float
        Mass for csp given the SFH input
    '''
    # Time array for sfh
    timek = np.logspace(-3, 2, 2000)
    # FIXME make this general for any sfh_func, no solution in mind yet
    sfh = sfh_func(timek, 10**(theta[3]), theta[5], theta[6], theta[7])
    ageval = 10**(theta[4])
    # Take only ages < age of galaxy for modeling
    sel = ages <= ageval
    # SFR interpolated from the finer sfh
    sfr = np.interp(ageval - ages, timek, sfh)
    # The weight is time between ages of each SSP
    weight = np.diff(np.hstack([0, ages])) * 1e9 * sfr
    # ages > ageval have zero weight
    weight[~sel] = 0
    # Cover the two cases where there is a fractional age range or not
    A = np.nonzero(ages <= ageval)[0][-1]
    B = np.nonzero(ages >= ageval)[0][0]
    if A == B:
        arr = np.dot(seds, weight)
        mass = np.sum(weight*masses)
    else:
        lw = ageval - ages[A]
        wei = lw * 1e9 * np.interp(ageval, timek, sfh)
        weight[B] = wei
        arr = np.dot(seds, weight)
        mass = np.sum(weight * masses)

    arr_dust = arr * np.exp(-1 * noll(wave, theta[0], theta[1], theta[2]))
    # Redshift
    csp = np.interp(wave, wave * (1+zobs), arr_dust)
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


def dust_absorption():
    pass


def dust_emmission():
    pass


def run_emcee():
    # Put SEDs for fitting into units at the given redshift
    pass


def plot_results():
    pass


def output_results():
    pass


def draw_uniform_dist(nsamples, start, end):
    return np.random.rand(nsamples)*(end-start) + start


def mock_data(args, ages, masses, wave, SSP, filter_matrix, filter_flag,
              nsamples=100):
    # Build fake theta, set z, mass, age to get sfh_a
    mock_list = ['mass', 'redshift', 'dust_tau', 'dust_delta', 'dust_bump',
                 'sfh_tau', 'sfh_a', 'sfh_b', 'sfh_c']
    theta = []
    for mock in mock_list:
        theta.append(draw_uniform_dist(nsamples,
                                       getattr(args, 'mock_' + mock)[0],
                                       getattr(args, 'mock_' + mock)[1]))
    zobs = draw_uniform_dist(nsamples, 1.9, 2.35)
    theta = np.array(theta).swapaxes(0, 1)

    y_model = []

    mass = []
    for thet in theta:
        spec, m = build_csp(thet, ages, SSP, masses, wave, zobs)
        mass.append(m)
        y = get_filter_fluxdensities(spec, filter_flag, filter_matrix)
        y.append()

    # calculate theta given input mass (need age or SFR)

    # get flux/mass back
    # output y, yerr, z, truth
    pass


def main(argv=None):
    # Get Inputs
    args = parse_args(argv)

    # Load Single Stellar Population model(s)
    ages, masses, wave, SSP = read_ssp(args)

    # Build Filter Matrix
    filter_matrix = build_filter_matrix(args, wave)

    # Make one instance of Mcsed for speed on initialization
    # Then replace the key variables each iteration for a given galaxy
    # Load sources for modelling
    if args.test:
        filter_flag = get_test_filters(args)
        y, yerr, z, flag, truth = mock_data(args, ages, masses, wave, SSP,
                                            filter_matrix, filter_flag)
    else:
        y, yerr, z, flag = read_input_file(args)

    
    for yi, ye, zi, fi in zip(y, yerr, z, flag):
        Mcsed()
    run_emcee()

    plot_results()
    
    output_results()
        
    
        
if __name__ == '__main__':
    main() 