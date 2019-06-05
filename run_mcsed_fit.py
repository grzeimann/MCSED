""" script for running MCSED

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

from __future__ import absolute_import
import sys
import argparse as ap
import numpy as np
import os.path as op
import logging
import config
from ssp import read_ssp
from astropy.io import fits
from astropy.table import Table, vstack
from mcsed import Mcsed
from distutils.dir_util import mkpath
from cosmology import Cosmology
# WPBWPB -- can remove this import
import time

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


def str2bool(v, log):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        log.warning('Could not interpret "metallicity" argument, by '
                    'default it will be set to False')
        return False


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

    # Avoid an infinite loop between parallel and series functions
    # only set to True if parallel mode already initiated 
    # from within run_mcsed_parallel.py script
    already_parallel = False
    if type(argv) is list:
        if '--already_parallel' in argv:
            argv.remove('--already_parallel')
            already_parallel = True

    # Pass "count" keyword -- only used in test mode in parallel
    count = 0
    if type(argv) is list:
        if '--count' in argv:
            indx = argv.index('--count')
            count = int( argv[indx+1] )
            del argv[indx: indx+2]

    parser = ap.ArgumentParser(description="MCSED",
                               formatter_class=ap.RawTextHelpFormatter)

    parser.add_argument("-f", "--filename",
                        help='''File to be read for galaxy data''',
                        type=str, default=None)

    parser.add_argument("-s", "--ssp",
                        help='''SSP Models, default fsps''',
                        type=str, default=None)

    parser.add_argument("-z", "--metallicity",
                        help='''Fixed metallicity for SSP models (0.02 is solar), False if free parameter''',
                        type=str, default=None)

    parser.add_argument("-i", "--isochrone",
                        help='''Isochrone for SSP model, e.g. padova''',
                        type=str, default=None)

    parser.add_argument("-sfh", "--sfh",
                        help='''Star formation history, e.g. constant''',
                        type=str, default=None)

    parser.add_argument("-dl", "--dust_law",
                        help='''Dust law, e.g. calzetti''',
                        type=str, default=None)

    parser.add_argument("-t", "--test",
                        help='''Test script with fake data''',
                        action="count", default=0)

    parser.add_argument("-tf", "--test_field",
                        help='''Test filters will match the given field''',
                        type=str, default='cosmos')

    parser.add_argument("-o", "--output_filename",
                        help='''Output filename for given run''',
                        type=str, default='test.dat')

    parser.add_argument("-nw", "--nwalkers",
                        help='''Number of walkers for EMCEE''',
                        type=int, default=None)

    parser.add_argument("-ns", "--nsteps",
                        help='''Number of steps for EMCEE''',
                        type=int, default=None)

    parser.add_argument("-an", "--add_nebular",
                        help='''Add Nebular Emission''',
                        action="count", default=0)

    parser.add_argument("-lu", "--logU",
                        help='''Ionization Parameter for nebular gas''',
                        type=float, default=None)

    parser.add_argument("-fd", "--fit_dust_em",
                        help='''Fit Dust Emission''',
                        action="count", default=0)

    parser.add_argument("-pf", "--phot_floor_error",
                        help='''Error floor for photometry''',
                        type=float, default=None)

    parser.add_argument("-ef", "--emline_floor_error",
                        help='''Error floor for emission lines''',
                        type=float, default=None)

    parser.add_argument("-no", "--nobjects",
                        help='''Number of test objects''',
                        type=int, default=None)

    parser.add_argument("-p", "--parallel",
                        help='''Running in parallel?''',
                        action="count", default=0)

    # Initialize arguments and log
    args = parser.parse_args(args=argv)
    args.log = setup_logging()

# WPBWPB delete
# add nebular becomes true no matter what argument is passed...
#    print(args.add_nebular)
#    print('args.test, getattr:  %s,  %s' % (args.test, getattr(args, 'test')))

    # Use config values if none are set in the input
    # Use command line arguments, if given. Else, refer to config file
#WPBWPB delete:
# i can combine this and reading from config into one list. list all desired variables
# and if cannot be passed on command line, it will already defer to config file
# unused from parser args... filename, output_filename, test, test_field, parallel, count
    arg_inputs = ['ssp', 'metallicity', 'isochrone', 'sfh', 'dust_law',
                  'nwalkers', 'nsteps', 'add_nebular', 'logU', 'fit_dust_em',
                  'phot_floor_error', 'emline_floor_error', 'nobjects',
                  'dust_em', 'Rv', 'EBV_stars_gas', 'wave_dust_em',
                  'emline_list_dict', 'emline_factor', 'use_input_data',
                  'metallicity_mass_relationship', 'metallicity_dict',
                  'filt_dict', 'catalog_filter_dict', 'filter_matrix_name',
                  'catalog_maglim_dict',
                  'output_dict', 'param_percentiles', 'reserved_cores']
    for arg_i in arg_inputs:
        try:
            if getattr(args, arg_i) in [None, 0]:
                setattr(args, arg_i, getattr(config, arg_i))
        except AttributeError:
            setattr(args, arg_i, getattr(config, arg_i))

    # If a test field is specified on the command line, initialize test mode
    if '-tf' in argv:
        args.test = True

#WPBWPB delete
#    print('test, fitdustem, nebular, parallel: %s, %s, %s, %s' % (args.test, args.fit_dust_em, args.add_nebular, args.parallel) )


#    args.fit_dust_em = str2bool(str(args.fit_dust_em), args.log)
#WPBWPB delete
#    print(args.fit_dust_em)
#    print(type(args.fit_dust_em))
#    print(args.add_nebular)

    # Set metallicity as free or fixed parameter
# WPBWPB: should 0,1 be boolean arguments on command line? boolean command line args only require flag, not flag+value -- metallicity is only ambiguous case...
    try:
        if args.metallicity not in ['0','1']:
            args.metallicity = float(args.metallicity)
    except ValueError:
        args.metallicity = str2bool(str(args.metallicity),args.log)
        print((args.metallicity, type(args.metallicity)))
        if args.metallicity:
            print("Fixing metallicity at z = 0.0077")
            args.metallicity = 0.0077

    # Avoid an infinite loop between parallel and series functions
    if already_parallel:
        args.already_parallel = True
    else:
        args.already_parallel = False

    # pass "count" keyword -- only used in test mode in parallel
    args.count = count

    # Determine whether emission lines should be carried
    # Collapses the emission line SSP grid to include only the desired lines
    # and speeds up the CSP construction
    if not args.add_nebular:
        args.use_emline_flux = False
    elif (args.emline_list_dict in [None, {}]):
        args.use_emline_flux = False
    elif (not args.test) & (not args.use_input_data):
        args.use_emline_flux = False
    else:
        args.use_emline_flux = True

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
# WPBWPB: either build from scratch each time, or check whether same dimensions
# only takes <0.5 sec to build -- might take longer just to load the thing
    if op.exists(args.filter_matrix_name):
        return np.loadtxt(args.filter_matrix_name)
    else:
# WPBWPB -- could remove
        start_time = time.time()
        nfilters = len(args.filt_dict)
        Fil_matrix = np.zeros((len(wave), nfilters))
        for i in np.arange(nfilters):
            wv, through = np.loadtxt(op.join('FILTERS', args.filt_dict[i]),
                                     unpack=True)
            new_through = np.interp(wave, wv, through, 0.0, 0.0)
            S = np.sum(new_through)
            if S == 0.:
                S = 1.
            Fil_matrix[:, i] = new_through / S
        np.savetxt(args.filter_matrix_name, Fil_matrix)
# WPBWPB -- could remove 
        ellapsed_time = time.time() - start_time
        print('Time to build filter matrix: %s sec' % ellapsed_time)
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
        if i in args.catalog_filter_dict[args.test_field]:
            filter_flag[i] = True
    return filter_flag


def get_maglim_filters(args):
    '''Used in test mode, this function loops through args.filt_dict and
    retrieves the 5-sigma magnitude limit for each filter (AB), and returns
    the appropriate microjansky 1-sigma error.

    Parameters
    ----------
    args : class
        The args class is carried from function to function with information
        from command line input and config.py

    Returns
    -------
    photerror : numpy array (float)
        Explained above.
    '''
    nfilters = len(args.filt_dict)
    photerror = np.zeros((nfilters,), dtype=float)
    for i in args.filt_dict.keys():
        if i in args.catalog_filter_dict[args.test_field]:
            maglim = args.catalog_maglim_dict[args.test_field][i]
            photerror[i] = 10**(-0.4 * (maglim - 23.9)) / 5.
    return photerror


def read_input_file(args):
    '''This function reads a very specific input file and joins it with
    archived 3dhst catalogs.  The input file should have the following columns:
    WPB
    FIELD, ID, Z

    Parameters
    ----------
    args : class
        The args class is carried from function to function with information
        from command line input and config.py

    Returns
    -------
    WPB: unit check: f_nu in microJy ?
    y : numpy array (2 dim)
        Photometric magnitudes from the 3DHST survey for each input source
    yerr : numpy array (2 dim)
        Photometric errors in magnitudes
    z : numpy array (1 dim)
        Redshift from the file returned as a numpy array
    flag : numpy array (2 dim)
        Flag set to True for filters in the catalog_filter_dict in config.py
    em : Astropy Table (2 dim)
        Emission line fluxes in ergs / cm2 / s
        Desired lines are read from dictionary in config.py
    emerr : Astropy Table (2 dim)
        Emission line errors in ergs / cm2 / s 
    '''
    # WPB: check if ID is of form skelton: if yes, grab from catalogs
    # else, check input units and convert appropriately
    F = Table.read(args.filename, format='ascii')
    nobj = len(F['field'])
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

    # keep track of which columns from the input file are utilized
    Fcols = list(F.colnames)
## WPBWPB delete
#    print('this is Fcols:   '+str(Fcols))

    # redshift array
    z = F['z']
    Fcols.remove('z')

    # convert from mag_zp = 25 to microjanskies (mag_zp = 23.9)
    fac = 10**(-0.4*(25.0-23.9))

    # read in emission line fluxes, if provided
    emline_fill_value = -99 # null value, should not be changed
    if args.use_emline_flux:
        em, emerr = Table(), Table()
        for emline in args.emline_list_dict.keys():
            colname, ecolname = '%s_FLUX' % emline, '%s_ERR' % emline
            if colname in Fcols:
                em_arr = np.array(F[colname]  * args.emline_factor)
                emerr_arr = np.max([F[ecolname], 
                                    args.emline_floor_error*F[colname]],0)
                emerr_arr *= args.emline_factor

                # account for objects with null measurements
                null_idx = np.where(np.array(F[colname])==emline_fill_value)[0]
                em_arr[null_idx] = emline_fill_value
                emerr_arr[null_idx] = emline_fill_value

                em[colname]     = em_arr 
                emerr[ecolname] = emerr_arr 

                Fcols = [c for c in Fcols if c not in [colname, ecolname]]
                print('Reading %s line fluxes from input file' % emline)
## WPBWPB delete
#                print('this is Fcols:   '+str(Fcols))
            else: 
## WPBWPB make a decision: all fill values, or remove from dictionary?
#                em[colname]     = np.full(len(F), emline_fill_value)
#                emerr[ecolname] = np.full(len(F), emline_fill_value)
                del args.emline_list_dict[emline]
    else:
        em    = np.full((len(F),2), emline_fill_value)
        emerr = np.full((len(F),2), emline_fill_value)

    for i, datum in enumerate(F):
        loc = datum[0].lower()
        for j, ind in enumerate(args.filt_dict.keys()):
            if ind in args.catalog_filter_dict[loc].keys():
                colname = "f_"+args.catalog_filter_dict[loc][ind]
                ecolname = "e_"+args.catalog_filter_dict[loc][ind]
            else:
                y[i, j] = 0.0
                yerr[i, j] = 0.0
                flag[i, j] = False
                continue
            if colname in field_dict[loc].columns.names:
                fi = field_dict[loc].data[colname][int(datum[1])-1]
                fie = field_dict[loc].data[ecolname][int(datum[1])-1]
                if (fi > -99.):
                    y[i, j] = fi*fac
                    # use a floor error if necessary
                    if fi != 0:
                        yerr[i, j] = np.abs(np.max([args.phot_floor_error,
                                            np.abs(fie/fi)]) * fi * fac)
                    else:
                        yerr[i, j] = 0.0
                        flag[i, j] = False
                    flag[i, j] = True
                else:
                    y[i, j] = 0.0
                    yerr[i, j] = 0.0
                    flag[i, j] = False
            else:
                y[i, j] = 0.0
                yerr[i, j] = 0.0
                flag[i, j] = False

    # warn of any unused columns from the input file
    if (Fcols!=[]) & (args.use_input_data):
        print('*CAUTION* unread columns in the input file:')
        print(Fcols) 

    # WPBWPB: adjust, clarify the column names
    return y, yerr, z, flag, F['obj_id'], F['field'], em, emerr


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


def mock_data(args, mcsed_model, nsamples=5, phot_error=0.05):
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
    WPB: unit check: f_nu in microJy ?
WPBWPB: modify, document all outputs
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
# WPB: modify, redshift range
    np.random.seed()
    thetas = mcsed_model.get_init_walker_values(num=nsamples, kind='ball')
    zobs = draw_uniform_dist(nsamples, 1.9, 2.35)
    params, y, yerr, true_y = [], [], [], []

    # add emission line fluxes
    emline_fill_value = -99 # null value, should not be changed
    if args.use_emline_flux: 
        em, emerr = Table(), Table()
        emlines = args.emline_list_dict.keys()
    else:
        em      = np.full( (nsamples,2), emline_fill_value)
        emerr   = np.full( (nsamples,2), emline_fill_value)
    for theta, z in zip(thetas, zobs):
        mcsed_model.set_class_parameters(theta)
        mcsed_model.set_new_redshift(z)
# WPBWPB delete
        print('this is the redshift: z')
        mcsed_model.spectrum, mass = mcsed_model.build_csp()
        # simulate emission line fluxes
        # true modeled fluxes are stored in mcsed_model.linefluxCSPdict
        if args.use_emline_flux: 
            em_loc, emerr_loc = Table(), Table()
            for emline in emlines:
## WPBWPB delete
#                print('emerr_loc: %s' % emerr_loc)
                colname, ecolname = '%s_FLUX' % emline, '%s_ERR' % emline
                model_lineflux = mcsed_model.linefluxCSPdict[emline]
## WPBWPB delete
#                print('err and type: %s, %s' % (model_lineflux * args.emline_floor_error, type(model_lineflux * args.emline_floor_error)))
#                print('colname, ecolname, and types: %s, %s, %s, %s' % (colname, ecolname, type(colname), type(ecolname)))
                emerr_loc[ecolname] = [model_lineflux * args.emline_floor_error]
## WPBWPB delete
#                perturb = float(emerr_loc[ecolname]*np.random.randn(1))
#                print(perturb) 
                em_loc[colname] = [model_lineflux 
                                + float(emerr_loc[ecolname]*np.random.randn(1))]
# WPBWPB delete
            print('these are the true and modeled line fluxes and errors:')
            print(mcsed_model.linefluxCSPdict.values())
            print(em_loc)
            print(emerr_loc)
            em = vstack([em, em_loc])
            emerr = vstack([emerr, emerr_loc])
## WPBWPB adjust log info.....
#        args.log.info('%0.2f, %0.2f' % (hlims[-1]*1e17, np.log10(mass)))
        f_nu = mcsed_model.get_filter_fluxdensities()
        if args.test_field in args.catalog_maglim_dict.keys():
            f_nu_e = get_maglim_filters(args)[mcsed_model.filter_flag]
            f_nu_e = np.max([f_nu_e, f_nu * phot_error], axis=0)
        else:
            f_nu_e = f_nu * phot_error
        y.append(f_nu_e*np.random.randn(len(f_nu)) + f_nu)
        yerr.append(f_nu_e)
        true_y.append(f_nu)
        params.append(list(theta) + [np.log10(mass)])

    return y, yerr, zobs, params, true_y, em, emerr


def main(argv=None, ssp_info=None):
    '''
    Execute the main functionality of MCSED

    Test mode: "python run_mcsed_fit.py -t"

    Live mode: "python run_mcsed_fit.py -f test_data.dat"

    For a "live" run, the key input ("-f") is a file with three columns:
    FIELD ID REDSHIFT

    The field options are: cosmos, goodsn, goodss, or aegis
    The id is the skelton photometric id for the given field
    The redshift is fixed in the fitting
    '''

    # Make output folder if it doesn't exist
    mkpath('output')

    # Get Inputs
    if argv == None:
        argv = sys.argv
        argv.remove('run_mcsed_fit.py')

    args = parse_args(argv)

## WPBWPB delete
#    print('this is argv from *fit.py.main():')
#    print(argv)
#    print(vars(args).keys())
#    print(args)
#    print(type(argv))
#    return


    # Run in parallel, if appropriate
    # (and if not already called from run_mcsed_parallel.py)
    if (args.parallel) & (not args.already_parallel):
## WPBWPB delete
#        print('just now starting parallel mode from within series function')
#        print('this is argv from *fit.py.main():')
#        print(argv)

        import run_mcsed_parallel
        run_mcsed_parallel.main_parallel(argv=argv)
        return   

##WPBWPB delete
#    print('running in series. here are args:')
#    print('this is argv from *fit.py.main():')
#    print(argv)
#    print(vars(args).keys())
#    print(args)
#    print(type(argv))

## WPBWPB delete
#    print('testing if boolean command line args default to True/False. this is add_nebular and fit_dust_emission:')
#    print((args.add_nebular, args.fit_dust_em))

    # Load Single Stellar Population model(s)
    if ssp_info is None:
## WPB delete: if carrying nebular separately from rest of spectrum...
#        ages, masses, wave, SSP, met, nebspec = read_ssp(args)
#        # WPBWPB delete
#        np.savez('spectra_separate',wave=wave,age=ages,spec=SSP,nspec=nebspec)
#        return
#
## WPBWPB delete: save emission line SSP grid
#        ages, masses, wave, SSP, met, nebSSP, linewave, lineSSP = read_ssp(args)
#        np.savez('emlineSSP', wave=linewave, age=ages, met=met, spec=lineSSP)
#        return
## WPB delete: if carrying them together, as originally done
## in that case, use the original (unmodified) ssp.py file
#        ages, masses, wave, SSP, met = read_ssp(args)
#        np.savez('spectra_together',wave=wave,age=ages,spec=SSP)
#        return

        args.log.info('Reading in SSP model')
        ages, masses, wave, SSP, met, linewave, lineSSP = read_ssp(args)

    else:
        ages, masses, wave, SSP, met, linewave, lineSSP = ssp_info

## WPBWPB delete
#    print((wave.shape, SSP.shape))
#    print(ages)
#    return


# WPBWPB: issue: want to read in the input file before building filter matrix,
# or build in another way -- need to know whether columns will be added
# also relevant for the emission line dictionary -- may exclude some lines

    # Build Filter Matrix
    filter_matrix = build_filter_matrix(args, wave)

## WPBWPB delete
#    print('filter_matrix.shape: %s' % filter_matrix.shape)
#    return

    # Make one instance of Mcsed for speed on initialization
    # Then replace the key variables each iteration for a given galaxy
    mcsed_model = Mcsed(filter_matrix, SSP, linewave, lineSSP, ages, 
                        masses, met, wave, args.sfh,
                        args.dust_law, args.dust_em, nwalkers=args.nwalkers,
                        nsteps=args.nsteps)

    # Communicate emission line measurement preferences
    mcsed_model.use_emline_flux = args.use_emline_flux
    mcsed_model.emline_dict = args.emline_list_dict

## WPBWPB delete
#    print('This is emline dict:')
#    print(args.emline_list_dict.keys())

## WPBWPB delete
#    return
#    # WPB delete -- modify to following section and uncomment
#    print(mcsed_model.dust_abs_class.Av)
#    print(mcsed_model.dust_abs_class.Rv)
#    print(mcsed_model.dust_abs_class.evaluate(np.array([5000])))
#    mcsed_model.dust_abs_class.Rv = args.Rv
#    print(mcsed_model.dust_abs_class.Rv)
#    print(mcsed_model.dust_abs_class.evaluate(np.array([5000])))

    # Adjust Rv in the dust attenuation model, if specified in config file
    # Otherwise, use the default value for the requested dust law
    if args.Rv >= 0:
        mcsed_model.dust_abs_class.Rv = args.Rv
    else:
        args.Rv = mcsed_model.dust_abs_class.Rv

    # Adjust the relative attenuation between stars and gas in the dust model
    # E(B-V)_stars = EBV_stars_gas * E(B-V)_gas
    if args.EBV_stars_gas >= 0:
        mcsed_model.dust_abs_class.EBV_stars_gas = args.EBV_stars_gas
    else:
        args.EBV_stars_gas = mcsed_model.dust_abs_class.EBV_stars_gas

#WPBWPB delete: this is where I adjust additional class params from config.py options

    # Special case of fixed metallicity
    if args.metallicity:
        mcsed_model.ssp_class.fix_met = True
        mcsed_model.ssp_class.met = args.metallicity
    else:
        mcsed_model.ssp_class.fix_met = False

    if not args.fit_dust_em:
        mcsed_model.dust_em_class.fixed = True
    else:
        mcsed_model.dust_em_class.fixed = False
# WPBWPB: maybe make more concise -- mcsed_model.dust_em_class.fixed = not args.fit_dust_em

    # Build names for parameters and labels for table
    names = mcsed_model.get_param_names()
    names.append('Log Mass')
    percentiles = args.param_percentiles 
    # WPB field/id
    labels = ['Field', 'ID', 'z']
    for name in names:
        labels = labels + [name + '_%02d' % per for per in percentiles]
    formats = {}
    for label in labels:
        formats[label] = '%0.3f'

    # If in test mode add truth values for table labels
    if args.test:
        for name in names:
            labels.append(name + '_truth')
            formats[labels[-1]] = '%0.3f'
    # WPB field/id
    formats['Field'], formats['ID'] = ('%s', '%05d')

##WPBWPB delete
#    print(formats)

    # WPB field/id
    mcsed_model.table = Table(names=labels, dtype=['S7', 'i4'] +
                              ['f8']*(len(labels)-2))

    # MAIN FUNCTIONALITY
    if args.test:
## WPBWPB delete
#        print('I"m in test mode!')
        fl = get_test_filters(args)
        mcsed_model.filter_flag = fl * True
        default = mcsed_model.get_params()
        y, yerr, z, truth, true_y, em, emerr = mock_data(args, mcsed_model,
                                                    phot_error=args.phot_floor_error,
                                                    nsamples=args.nobjects)

        cnts = np.arange(args.count, args.count + len(z))
# WPB: change emline input (?)
        for yi, ye, zi, tr, ty, cnt, emi, emie in zip(y, yerr, z, truth, true_y, 
                                               cnts, em, emerr):
            mcsed_model.input_params = tr
            mcsed_model.filter_flag = fl * True
            mcsed_model.set_class_parameters(default)
            mcsed_model.data_fnu = yi
            mcsed_model.data_fnu_e = ye
            mcsed_model.true_fnu = ty
            mcsed_model.set_new_redshift(zi)
# WPBWPB: modify? true or perturbed emission line flux?
            mcsed_model.data_emline = emi
            mcsed_model.data_emline_e = emie

            # Remove filters containing Lyman-alpha (and those blueward)
            mcsed_model.remove_waverange_filters(0., 1216., restframe=True)
            # Remove filters dominated by dust emission, if applicable
            if not args.fit_dust_em:
                mcsed_model.remove_waverange_filters(args.wave_dust_em*1e4,1e10,
                                                     restframe=True)

            mcsed_model.fit_model()

            if args.output_dict['sample plot']:
                mcsed_model.sample_plot('output/sample_fake_%05d_%s' % (cnt, let))
            if args.output_dict['triangle plot']:
                mcsed_model.triangle_plot('output/triangle_fake_%05d_%s_%s_%s' % (cnt, let, args.sfh, args.dust_law))

            mcsed_model.table.add_row(['Test', cnt, zi] + [0.]*(len(labels)-3))
            last = mcsed_model.add_fitinfo_to_table(percentiles)
            mcsed_model.add_truth_to_table(tr, last)
            print(mcsed_model.table)
    else:
    # WPB field/id
        y, yerr, z, flag, objid, field, em, emerr = read_input_file(args)
##WPBWPB delete
#        print(read_input_file(args))
#        print(em)
## WPBWPB delete
#        print('This is emline dict:')
#        print(args.emline_list_dict.keys())
#        return

        # Update emission line dictionary (remove lines absent in input file)
        mcsed_model.emline_dict = args.emline_list_dict

        iv = mcsed_model.get_params()
# WPBWPB delete
#        print(iv)
        #return
        for yi, ye, zi, fl, oi, fd, emi, emie in zip(y, yerr, z, flag, objid,
                                                   field, em, emerr):
            mcsed_model.filter_flag = fl
            mcsed_model.set_class_parameters(iv)
# WPBWPB delete
#            return
            mcsed_model.data_fnu = yi[fl]
            mcsed_model.data_fnu_e = ye[fl]
            mcsed_model.set_new_redshift(zi)
            mcsed_model.data_emline = emi
            mcsed_model.data_emline_e = emie

## WPB delete
#            fwave = mcsed_model.get_filter_wavelengths()
#            print('these are filter wavelengths:')
#            print(np.sort(fwave))
#            dummy = mcsed_model.remove_waverange_filters(79560.,1e9,restframe=False)
#            fwave = mcsed_model.get_filter_wavelengths()
#            print(np.sort(fwave))
#            return

            # Remove filters containing Lyman-alpha (and those blueward)
            mcsed_model.remove_waverange_filters(0., 1216., restframe=True)
            # Remove filters dominated by dust emission, if applicable
            if not args.fit_dust_em:
                mcsed_model.remove_waverange_filters(args.wave_dust_em*1e4,1e10, 
                                                     restframe=True)

## WPB delete
#            fwave = mcsed_model.get_filter_wavelengths()
#            print('these are filter wavelengths:')
#            print(np.sort(fwave))
#            return

            mcsed_model.fit_model()
##WPBWPB delete
#            print('I"ve reached this point')
    # WPB field/id
            if args.output_dict['sample plot']:
                mcsed_model.sample_plot('output/sample_%s_%05d' % (fd, oi), 
                                        imgtype = args.output_dict['image format'])
            if args.output_dict['triangle plot']:
                mcsed_model.triangle_plot('output/triangle_%s_%05d_%s_%s' %
                                          (fd, oi, args.sfh, args.dust_law),
                                          imgtype = args.output_dict['image format'])
            mcsed_model.table.add_row([fd, oi, zi] + [0.]*(len(labels)-3))
            names = mcsed_model.get_param_names()
            names.append('Log Mass')
            names.append('Ln Prob')
            if args.output_dict['fitposterior']:
                T = Table(mcsed_model.samples, names=names)
                T.write('output/fitposterior_%s_%05d_%s.dat' % (fd, oi, args.sfh),
                        overwrite=True, format='ascii.fixed_width_two_line')
            if args.output_dict['bestfitspec']:
                T = Table([mcsed_model.wave, mcsed_model.medianspec],
                          names=['wavelength', 'spectrum'])
                T.write('output/bestfitspec_%s_%05d_%s.dat' % (fd, oi, args.sfh),
                        overwrite=True, format='ascii.fixed_width_two_line')
            if args.output_dict['bestfitflux']:
                T = Table([mcsed_model.fluxwv, mcsed_model.fluxfn],
                          names=['wavelength', 'fluxdensity'])
                T.write('output/bestfitflux_%s_%05d_%s.dat' % (fd, oi, args.sfh),
                        overwrite=True, format='ascii.fixed_width_two_line')
            if args.output_dict['observedflux']:
                T = Table([mcsed_model.fluxwv, mcsed_model.data_fnu,
                           mcsed_model.data_fnu_e],
                           names=['wavelength', 'fluxdensity', 'fluxdensityerror'])
                T.write('output/observedflux_%s_%05d_%s.dat' % (fd, oi, args.sfh),
                        overwrite=True, format='ascii.fixed_width_two_line')
            last = mcsed_model.add_fitinfo_to_table(percentiles)
            print(mcsed_model.table)
    if args.parallel:
        return [mcsed_model.table, formats]
    else:
        if args.output_dict['parameters']:
            mcsed_model.table.write('output/%s' % args.output_filename,
                                    format='ascii.fixed_width_two_line',
                                    formats=formats, overwrite=True)
        if args.output_dict['settings']:
            filename = open('output/%s.args' % args.output_filename, 'w')
            filename.write( str( vars(args) ) )
            filename.close()
if __name__ == '__main__':
    main()


