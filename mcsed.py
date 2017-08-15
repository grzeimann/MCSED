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
import emcee
import config
from dust_abs import calzetti, noll
#from dust_em import draine_li
from sfh import double_powerlaw
from ssp import read_ssp

def parse_args(argv=None):
    # Arguments to parse include date range, inspection attribute, instrument,
    #           instrument_side
    parser = ap.ArgumentParser(description="MCSED",
                            formatter_class=ap.RawTextHelpFormatter)
                            
    parser.add_argument("-s","--ssp", 
                        help='''SSP Models, default fsps''',
                        type=str, default=None)

    parser.add_argument("-z","--metallicity", 
                        help='''Metallicity of the SSP models, 0.02 is solar''',
                        type=float, default=None)

    parser.add_argument("-i","--isochrone", 
                        help='''Isochrone for SSP model, e.g. padova''',
                        type=str, default=None)
                   
    args = parser.parse_args(args=argv)
    
    arg_inputs = ['ssp', 'metallicity']
    for arg_i in arg_inputs:
        if getattr(args, arg_i) is None:
            setattr(args, arg_i, getattr(config, arg_i))
    
    config_copy_list = ['metallicity_dict']
    for con_copy in config_copy_list:
        setattr(args, con_copy, getattr(config, con_copy))
        
    return args

def read_input_file():
    pass

def build_filter_matrix():
    pass

def get_filter_fluxdensities():
    pass
 
def build_csp():
    pass

def lnlike():
    pass

def testbounds(theta):
    return False

def lnprior(theta):
    flag = testbounds(theta)
    if not flag:
        return -np.inf
    else:
        return 0.0

def lnprob():
    pass

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

def main(argv=None):
    # Get Inputs
    args = parse_args(argv)
    
    # Load Single Stellar Population model(s)
    read_ssp(args)
    
        
if __name__ == '__main__':
    main() 