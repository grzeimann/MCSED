""" MCSED - ssp.py

Single Stellar Population module for loading models

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import os.path as op
import numpy as np
import sys


def read_ssp(args):
    ''' Read in SPS model and return ages, masses, wavelength, and spectra

    Parameters
    ----------
    args : class
        The args class from mcsed.parse_args()

    '''
    if args.ssp.lower() == 'fsps':
        return read_fsps(args)


def read_fsps(args):
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
        wavelength for each spectrum
    spec : numpy array (2 dim)
        Spectra in f_nu (ergs/s/cm^2/Hz) at 10pc
    '''
    pc10 = 10. * 3.08567758e18
    solar_microjansky = 3.826e33 * 1e29 / (4. * np.pi * pc10**2)
    filename = op.join('SSP', 'fsps_%s_%0.4f.spec' % (args.isochrone,
                                                      args.metallicity))
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
                    spec.append(np.array(l, dtype=float) * solar_microjansky)
            cnt += 1
    sel = np.where((wave > 500) * (wave < 1e5))[0]
    spec = np.array(spec).swapaxes(0, 1)[sel, :]
    wave = wave[sel]
    nsel = np.where(np.abs(np.diff(wave)-0.9) < 0.5)[0]
    ndw = 8.
    nw = np.arange(wave[nsel[0]], wave[nsel[-1]+1], ndw)
    nwb = np.hstack([nw, wave[nsel[-1]+1]])
    nwave = np.hstack([wave[:nsel[0]], nw, wave[(nsel[-1]+1):]])
    nspec = np.zeros((len(nwave), spec.shape[1]))
    for i, sp in enumerate(spec.swapaxes(0, 1)):
        hsp = (np.histogram(wave[nsel], nwb, weights=sp[nsel])[0] /
               np.histogram(wave[nsel], nwb)[0])
        nspec[:, i] = np.hstack([sp[:nsel[0]], hsp, sp[(nsel[-1]+1):]])
    return (10**(np.array(ages)-9.), 10**np.array(masses), nwave, nspec)
