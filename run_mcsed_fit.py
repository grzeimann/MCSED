""" Script for running MCSED

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import argparse as ap
import numpy as np
import os.path as op
import logging
import config
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
    arg_inputs = ['ssp', 'metallicity', 'isochrone']
    for arg_i in arg_inputs:
        if getattr(args, arg_i) is None:
            setattr(args, arg_i, getattr(config, arg_i))

    # Copy list of config values to the args class
    config_copy_list = ['metallicity_dict', 'filt_dict', 'catalog_filter_dict',
                        'test_filter_dict', 'filter_matrix_name', 'sfh',
                        'dust_law']

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
            wv, through = np.loadtxt(op.join('FILTERS', args.filt_dict[i]),
                                     unpack=True)
            new_through = np.interp(wave, wv, through, 0.0, 0.0)
            Fil_matrix[:, i] = new_through / np.sum(new_through)
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


def draw_uniform_dist(nsamples, start, end):
    ''' Draw random samples from a uniform distribution

    Parameters
    ----------
    nsamples : int
        Number of draws
    start : float
        lower bound
    end : float
        higher bound

    Returns
    -------
    uniform_sample : numpy array (1 dim)
        randomly drawn variables from a uniform distribution
    '''
    return np.random.rand(nsamples)*(end-start) + start


def draw_gaussian_dist(nsamples, means, sigmas):
    ''' Draw random samples from a normal distribution

    Parameters
    ----------
    nsamples : int
        Number of draws
    means : numpy array or list
        Average values to draw from
    sigmas : numpy array or list
        Standard deviation of the return distributions

    Returns
    -------
    normal_sample : numpy array (1 dim)
        randomly drawn variables from a normal distribution
    '''
    m = len(means)
    N = np.random.randn(nsamples * m).reshape(nsamples, m)
    return sigmas * N + means


def mock_data(args, mcsed_model, nsamples=5, phot_error=0.1):
    ''' Create mock data to test quality of MCSED fits

    Parameters
    ----------
    args : class
        The args class is carried from function to function with information
        from command line input and config.py
    mcsed_model : class
        Mcsed class for building fake galaxies given input thetas

    Returns
    -------
    y : numpy array (2 dim)
        Photometric magnitudes for mock galaxies
    yerr : numpy array (2 dim)
        Photometric errors in magnitudes
    z : numpy array (1 dim)
        Redshift for mock galaxies
    truth : numpy array (2 dim)
        Mock input parameters for each fake galaxy, e.g. dust, sfh, mass
    '''
    # Build fake theta, set z, mass, age to get sfh_a
    thetas = mcsed_model.get_init_walker_values(num=nsamples)
    zobs = draw_uniform_dist(nsamples, 1.9, 2.35)
    params, y, yerr = [], [], []
    for theta, z in zip(thetas, zobs):
        mcsed_model.set_class_parameters(theta)
        mcsed_model.set_new_redshift(z)
        mcsed_model.spectrum, mass = mcsed_model.build_csp()
        f_nu = mcsed_model.get_filter_fluxdensities()
        f_nu_e = f_nu * phot_error
        y.append(f_nu)
        yerr.append(f_nu_e)
        params.append(list(theta) + [np.log10(mass)])
        #wv = mcsed_model.get_filter_wavelengths()
        #import matplotlib.pyplot as plt
        #plt.plot(mcsed_model.wave, mcsed_model.spectrum)
        #plt.scatter(wv, f_nu, color='r')
        #plt.axis([4000,30000,0.01,100])
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.savefig('test.png')
        #import sys
        #sys.exit(1)
    return y, yerr, zobs, params


def main(argv=None):
    # Get Inputs
    args = parse_args(argv)

    # Load Single Stellar Population model(s)
    ages, masses, wave, SSP = read_ssp(args)
    # Build Filter Matrix
    filter_matrix = build_filter_matrix(args, wave)

    # Make one instance of Mcsed for speed on initialization
    # Then replace the key variables each iteration for a given galaxy
    mcsed_model = Mcsed(filter_matrix, SSP, ages, masses, wave, args.sfh,
                        args.dust_law)

    if args.test:
        mcsed_model.filter_flag = get_test_filters(args)
        y, yerr, z, truth = mock_data(args, mcsed_model)
        cnt = 0
        for yi, ye, zi, tr in zip(y, yerr, z, truth):
            mcsed_model.input_params = tr
            mcsed_model.set_class_parameters(tr[:-1])
            mcsed_model.data_fnu = yi
            mcsed_model.data_fnu_e = ye
            mcsed_model.set_new_redshift(zi)
            mcsed_model.fit_model()
            mcsed_model.triangle_plot('output/triangle_fake_%04d' % cnt)
            cnt += 1
    else:
        y, yerr, z, flag = read_input_file(args)


if __name__ == '__main__':
    main()
