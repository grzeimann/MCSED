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
import emcee
import config
from dust_abs import calzetti, noll
#from dust_em import draine_li
from sfh import double_powerlaw

def parse_args(argv=None):
    # Arguments to parse include date range, inspection attribute, instrument,
    #           instrument_side
    parser = ap.ArgumentParser(description="MCSED",
                            formatter_class=ap.RawTextHelpFormatter)
                            
    parser.add_argument("-s","--ssp", 
                        help='''SSP Models, default fsps''',
                        type=str, default='fsps')
                        
    args = parser.parse_args(args=argv)

    return args

def read_input_file():
    pass

def load_ssp():
    pass

def build_filter_matrix():
    pass

def get_filter_fluxdensities():
    pass
 
def build_csp():
    pass

def lnlike():
    pass

def lnprior():
    pass

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

def main():
    pass
    
if __name__ == '__main__':
    main() 