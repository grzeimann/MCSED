""" MCSED - dust.py


1) Dust Laws 
    a) Calzetti:
    b) Noll:

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np

def calzetti(wave):
    invwv = 1/(wave/1e4)
    sel1 = np.nonzero(wave<0.63e4)[0]
    sel2 = np.nonzero(np.logical_and(wave>=0.63e4,wave<2.2e4))[0]
    k1 = np.zeros(sel1.shape)
    k2 = np.zeros(sel2.shape)
    k1 = (2.659*(-2.156 + 1.509*invwv[sel1] -0.198*invwv[sel1]**2 
                 + 0.011*invwv[sel1]**3) + 4.05)
    k2 = 2.659*(-1.857 + 1.040*invwv[sel2]) + 4.05
    k = np.zeros(wave.shape)
    k[sel1] = k1
    k[sel2] = k2
    return k

def noll(wave, tau=0.0, delta=0.0, Eb=0.0):
    dellam = 350.
    lam0   = 2175.
    Dlam   = Eb * (wave*dellam)**2 / ((wave**2-lam0**2)**2+(wave*dellam)**2)
    taulam = tau / 4.05 * (calzetti(wave)+Dlam)*(wave/5500)**(delta)
    return taulam