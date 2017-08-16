""" MCSED - ssp.py


1) Single Stellar Populations

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import os.path as op
import numpy as np

def read_ssp(args):
    if args.ssp.lower() == 'fsps':
        return read_fsps(args)
    

def read_fsps(args):
    args.metallicity_dict[args.isochrone], args.metallicity
    filename = op.join('SPS', 'fsps_%s_%0.4f.spec' %(args.isochrone, met))
    cnt = 0
    ages = []
    masses = []
    spec = []
    with open(filename) as f:
        for lines in f:
            if cnt==9:
                wave = np.array(lines.split(' '))
            if cnt>9:
                l = lines.split(' ')
                if len(l) == 4:
                    ages.append(l[0])
                    masses.append(l[1])
                else:
                    spec.append(np.array(l))
    return np.array(ages), np.array(masses), wave, np.array(spec)