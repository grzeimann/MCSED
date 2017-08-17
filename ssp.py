""" MCSED - ssp.py


1) Single Stellar Populations

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import os.path as op
import numpy as np
import sys

def read_ssp(args):
    if args.ssp.lower() == 'fsps':
        return read_fsps(args)
    

def read_fsps(args):
    filename = op.join('SPS', 'fsps_%s_%0.4f.spec' %(args.isochrone, 
                                                     args.metallicity))
    if not op.exists(filename):
        print('Metallicity entered, %0.4f, does not match any for %s'
              'isochrones of the %s models' %(args.isochrone, args.ssp))
        print('Metallicity options [')
        for met in args.metallicity_dict[args.isochrone]:
            print('%0.4f ')
        print(']')
        sys.exit(1)
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