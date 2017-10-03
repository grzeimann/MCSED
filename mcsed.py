""" MCSED


1) CURRENT LIMITATIONS:
       A) Constant metallicity
       B) Dust Emission is AD HoC from Draine and Li (2007)
       
   OPTIONAL FITTED PARAMETERS:
       A) SFH
           a) For example: tau_sfh, age, a, b, c
       B) Dust law
           b) For example: tau_dust, delta, Eb
   
   DERIVED PARAMETERS:
       A) Stellar Mass
       B) Dust Emission

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import argparse as ap
import numpy as np
import os.path as op
import emcee
import config
from dust_abs import noll
#from dust_em import draine_li
from sfh import double_powerlaw
from ssp import read_ssp
from astropy.io import fits


def parse_args(argv=None):
    # Arguments to parse include ssp code, metallicity, isochrone choice, 
    #   whether this is real data or mock data
    parser = ap.ArgumentParser(description="MCSED",
                            formatter_class=ap.RawTextHelpFormatter)

    parser.add_argument("-f","--filename", 
                        help='''File to be read for galaxy data''',
                        type=str, default=None)
                            
    parser.add_argument("-s","--ssp", 
                        help='''SSP Models, default fsps''',
                        type=str, default=None)

    parser.add_argument("-z","--metallicity", 
                        help='''Metallicity of the SSP models, 0.02 is solar''',
                        type=float, default=None)

    parser.add_argument("-i","--isochrone", 
                        help='''Isochrone for SSP model, e.g. padova''',
                        type=str, default=None)

    parser.add_argument("-t","--test", 
                        help='''Test script with fake data''',
                        action="count", default=0)
                   
    args = parser.parse_args(args=argv)
    
    arg_inputs = ['ssp', 'metallicity']
    for arg_i in arg_inputs:
        if getattr(args, arg_i) is None:
            setattr(args, arg_i, getattr(config, arg_i))
    
    config_copy_list = ['metallicity_dict', 'mock_masses', 'mock_redshift',
                        'mock_dust_tau', 'mock_dust_delta', 'mock_dust_bump', 
                        'mock_sfh_tau', 'mock_sfh_b', 'mock_sfh_c',
                        'filt_dict', 'catalog_filt_dict',
                        'filter_matrix_name']
                        
    for con_copy in config_copy_list:
        setattr(args, con_copy, getattr(config, con_copy))
        
    return args

def build_filter_matrix(args, wave):
    if op.exists(args.filter_matrix_name):
        return np.loadtxt(args.filter_matrix_name)
    else:
        nfilters=len(args.filt_dict)
        Fil_matrix = np.zeros((len(wave),nfilters))
        for i in np.arange(nfilters):
            wv, through = np.loadtxt('FILTERS',args.filt_dict[i], unpack=True)
            new_through = np.interp(wave,wv,through,0.0,0.0)
            Fil_matrix[:,i] = new_through/np.sum(new_through)
        np.savetxt(args.filter_matrix_name, Fil_matrix)
        return Fil_matrix

def get_filter_flag(args):
    nfilters = len(args.filt_dict)
    filter_flag = np.ones((nfilters,))>0
    for i in args.filt_dict.keys():    
        try:
            args.catalog_filter_dict[i]
        except:
            filter_flag[i]=False
    return filter_flag   

def read_input_file(args):
    # read file
    F = np.loadtxt(args.filename)
    nobj = len(F)
    fields = ['aegis','cosmos','goodsn','goodss','uds']
    name_base = '_3dhst.v4.1.cat.FITS'
    field_dict = {}
    for field in fields:
        field_dict[field] = fits.open(op.join('3dhst_catalogs', 
                                              field+name_base))[1]
    nfilters = len(args.filt_dict)
    y = np.zeros((nobj,nfilters))
    yerr = np.zeros((nobj,nfilters))
    flag = np.zeros((nobj,nfilters),dtype=bool)
    z = F[:,2]
    for i,datum in enumerate(F):
        loc = datum[0].lower()
        for j,ind in enumerate(args.filt_dict.keys()):
            colname = "f_"+args.catalog_filter_dict[ind]
            ecolname = "e_"+args.catalog_filter_dict[ind]
            if colname in field_dict[loc].columns.names:
                y[i,j] = field_dict[loc].data[colname][int(datum[1])]
                yerr[i,j] = field_dict[loc].data[ecolname][int(datum[1])]
                flag[i,j] = True
            else:
                y[i,j] = 0.0
                yerr[i,j] = 0.0
                flag[i,j] = False
    return y, yerr, z, flag
            
    
def get_filter_fluxdensities(spec, filter_flag, filter_matrix):
    a = np.dot(spec, filter_matrix)
    return a[filter_flag]
 
def build_csp(theta, ages, seds, masses, wave, zobs):
    
    pass

def testbounds(theta):
    return False

def lnprior(theta):
    flag = testbounds(theta)
    if not flag:
        return -np.inf
    else:
        return 0.0

def lnlike(theta, ages, seds, masses, wave, zobs, y, yerr, filters, sigma_m):
    flag = testbounds(theta)
    if flag:
        return -np.inf, -np.inf
    else:
        flux, mass = build_csp(theta, ages, seds, masses, wave, zobs)
        model_y = get_filter_fluxdensities(wave, flux, filters)
        inv_sigma2     = 1.0/(yerr**2+(model_y*sigma_m)**2)
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
    pass

def plot_results():
    pass

def output_results():
    pass

def draw_uniform_dist(nsamples, start, end):
    return np.random.rand(nsamples)*(end-start) + start
    
def mock_data(args, nsamples=100):
    mock_list = ['mass', 'redshift', 'dust_tau', 'dust_delta', 'dust_bump', 
                 'sfh_tau', 'sfh_a', 'sfh_b', 'sfh_c']
    theta = []
    for mock in mock_list:
        theta.append(draw_uniform_dist(nsamples, getattr(args,'mock_'+mock)[0], 
                                       getattr(args,'mock_'+mock)[1]))
    # calculate theta given input mass (need age or SFR)
    '''
    tim       = ages / 1e9
    timek     = np.logspace(-3,2,2000)   
    sfh       = doublepowerlaw_sfh(timek,10**(theta[3]),theta[5],theta[6],theta[7])
    ageval    = 10**(theta[4])
    sel = tim <= ageval
    sfr       = np.interp(ageval - tim,timek,sfh)
    weight = np.diff(np.hstack([0,tim]))*1e9*sfr 
    weight[~sel]=0
    '''                                    
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
    filter_flag = get_filter_flag(args)
    
    # Load sources for modelling
    if args.test:
        y, yerr, z, flag, truth = mock_data(ages, masses, wave, SSP, 
                                            filter_matrix,
                                            filter_flag)
    else:
        y, yerr, z, flag = read_input_file(args)
        
    run_emcee()
    
    plot_results()
    
    output_results()
        
    
        
if __name__ == '__main__':
    main() 