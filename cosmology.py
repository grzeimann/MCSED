""" Cosmology Calculator

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np


class Cosmology:
    def __init__(self, omega_m=0.31, omega_l=0.69, omega_k=0.00, h0=0.69):
        ''' Initialize the class

        Inputs
        ------
        omega_m : float
            Matter density in the universe today
        omega_l : float
            Dark Energy density in the universe today
        omega_k : float
            Curvature density in the universe today
        h0 : float
            Hubble Constant in the universe today
        '''
        self.omega_m = omega_m
        self.omega_l = omega_l
        self.omega_k = omega_k
        self.h0 = h0
        self.c = 2.99792e18

    def luminosity_distance(self, z, stepsize=0.001):
        ''' Calculate the luminosity distance (cm)

        Inputs
        ------
        z : float
            Redshift
        stepsize : float
            Integration stepsize

        Returns
        -------
        d : float
            Luminosity distance (cm)
        '''
        zi = np.arange(0., z + stepsize, stepsize)
        E = np.sqrt(self.omega_m * (1 + zi)**3 + self.omega_k * (1 + zi)**2 +
                    self.omega_l)
        d = (1 + z) * self.c / self.h0 * np.sum(1. / E * stepsize)
        return d
