""" Cosmology Calculator

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np


class Cosmology:
    def __init__(self, omega_m=0.31, omega_l=0.69, omega_k=0.00, H_0=69):
        ''' Initialize the class

        Parameters
        ----------
        omega_m : float
            Matter density in the universe today
        omega_l : float
            Dark Energy density in the universe today
        omega_k : float
            Curvature density in the universe today
        H_0 : float
            Hubble Constant in the universe today
        '''
        self.omega_m = omega_m
        self.omega_l = omega_l
        self.omega_k = omega_k
        self.H_0 = H_0  # km / s / (10 pc)
        self.t_h = 9.78 * (100. / self. H_0) # Gyr
        self.c = 2.99792e5  # km / s

    def luminosity_distance(self, z, stepsize=0.001):
        ''' Calculate the luminosity distance

        Parameters
        ----------
        z : float
            Redshift
        stepsize : float
            Integration stepsize

        Returns
        -------
        d : float
            Luminosity distance (units of 10 pc)
        '''
        zi = np.arange(0., z, stepsize)
        E = np.sqrt(self.omega_m * (1 + zi)**3 + self.omega_k * (1 + zi)**2 +
                    self.omega_l)
        d = (1 + z) * self.c / self.H_0 * np.sum(1. / E * stepsize)
        return d * 1e5

    def lookback_time(self, z, stepsize=0.001):
        ''' Calculate the lookback time

        Parameters
        ----------
        z : float
            Redshift
        stepsize : float
            Integration stepsize

        Returns
        -------
        t : float
            Lookback time (units of Gyr)
        '''
        zi = np.arange(0., z, stepsize)
        E = np.sqrt(self.omega_m * (1 + zi)**3 + self.omega_k * (1 + zi)**2 +
                    self.omega_l)
        t = self.t_h * np.sum(stepsize / ((1. + zi) * E))
        return t
