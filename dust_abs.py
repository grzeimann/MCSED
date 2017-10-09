""" MCSED - dust.py


1) Dust Laws
    a) Calzetti:
    b) Noll:

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np


def calzetti(wave):
    ''' Calzetti et al. (2000) dust attenuation curve, k(wave)

    Input
    -----
    wave : numpy array (1 dim)
        wavelength

    Returns
    -------
    k : numpy array (1 dim)
        A(wave) = R_V * k(wave)
    '''
    invwv = 1/(wave/1e4)
    sel1 = np.nonzero(wave < 0.63e4)[0]
    sel2 = np.nonzero(np.logical_and(wave >= 0.63e4, wave < 2.2e4))[0]
    k1 = np.zeros(sel1.shape)
    k2 = np.zeros(sel2.shape)
    k1 = (2.659 * (-2.156 + 1.509 * invwv[sel1] - 0.198 * invwv[sel1]**2 +
          0.011 * invwv[sel1]**3) + 4.05)
    k2 = 2.659*(-1.857 + 1.040*invwv[sel2]) + 4.05
    k = np.zeros(wave.shape)
    k[sel1] = k1
    k[sel2] = k2
    return k


class noll:
    ''' Prescription for dust law comes from Noll et al. (2009), with constants
    defined in Kriek & Conroy (2013).  This dust attenuation law includes a
    bump at 2175A and a modified Calzetti et al. (2000) attenuation curve.

        A(wave) = frac{A_V}{4.05} (k'(wave) + D(wave))
                     left(frac{wave}{5500}right) ^delta
        D(wave) = frac{E_b (wave,dellam)^2 }{(wave^2-lam0^2)^2
                     + (wave,dellam)^2}
    '''
    def __init__(self, tau=0.0, delta=0.0, Eb=0.0):
        ''' Initialize Class

        Input
        -----
        tau : float
            Effective depth, e.g., Observed = True * exp**(-tau/4.05 * k(wave))
        delta : float
            Power for powerlaw modification of Calzetti curve
        Eb : float
            Strength of 2175A bump.  See equation above for the Drude profile,
            D(wave)
        '''
        self.tau = tau
        self.delta = delta
        self.Eb = Eb
        self.nparams = 3

    def set_parameters_from_list(self, input_list, start_value):
        ''' Set parameters from a list and a start_value

        Input
        -----
        input_list : list
            list of input parameters (could be much larger than number of
            parameters to be set)
        start_value : int
            initial index from list to read out parameters
        '''
        self.tau = input_list[start_value]
        self.delta = input_list[start_value+1]
        self.Eb = input_list[start_value+2]

    def evaluate(self, wave):
        ''' Evaluate Dust Law

        Input
        -----
        wave : numpy array (1 dim)
            wavelength

        Returns
        -------
        taulam : numpy array (1 dim)
            Effective optical depth as a function of wavelength
        '''
        dellam = 350.
        lam0 = 2175.
        Dlam = (self.Eb * (wave*dellam)**2 /
                          ((wave**2-lam0**2)**2+(wave*dellam)**2))
        taulam = (self.tau / 4.05 *
                  (calzetti(wave)+Dlam)*(wave/5500)**(self.delta))
        return taulam
