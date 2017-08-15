""" MCSED - ssp.py


1) Single Stellar Populations

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import os.path as op

def read_ssp(args):
    if args.ssp.lower() == 'fsps':
        return read_fsps(args)
    

def read_fsps(args):
    args.metallicity_dict[args.isochrone], args.metallicity
    filename = op.join('SPS', 'fsps_%s_%0.4f.spec' %(args.isochrone, met)
    return filename