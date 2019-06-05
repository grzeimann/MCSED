# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 11:09:49 2018

@author: gregz
"""

import numpy as np
from scipy.interpolate import LinearNDInterpolator


class DL07:
    ''' Prescription for dust emission comes from Draine & Li (2007) '''
    def __init__(self, umin=2.0, gamma=0.05, qpah=2.5, umin_lims=[0.1, 25.0],
                 gamma_lims=[0, 1.], qpah_lims=[0.47, 4.58], umin_delta=0.4,
                 gamma_delta=0.02, qpah_delta=1.0, fixed=True):
        ''' Initialize Class

        Parameters
        -----
        umin : float
            Minimum in distribution of incident intensity
        gamma : float
            Fraction of distribution in incident intensity greater than umin
        qpah : float
            Percentage of PAH emission
        '''
        self.umin = umin
        self.gamma = gamma
        self.qpah = qpah
        self.umin_lims = umin_lims
        self.gamma_lims = gamma_lims
        self.qpah_lims = qpah_lims
        self.umin_delta = umin_delta
        self.gamma_delta = gamma_delta
        self.qpah_delta = qpah_delta
        self.fixed = fixed
        self.get_dust_emission_matrix()

    def get_dust_emission_matrix(self):
        DE1, DE2 = (np.zeros((7, 22, 1001)), np.zeros((7, 22, 1001)))
        self.qpaharray = np.array([0.47, 1.12, 1.77, 2.50, 3.19, 3.90, 4.58])
        self.uminarray = np.array([0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8,
                                   1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0,
                                   8.0, 12.0, 15.0, 20.0, 25.0])
        Xgrid, Ygrid = np.meshgrid(self.qpaharray, self.uminarray)
        X = np.vstack([Xgrid.ravel(), Ygrid.ravel()]).swapaxes(0, 1)

        for j, i in enumerate(np.arange(0, 70, 10)):
            de = np.loadtxt('DUSTEMISSION/DL07/DL07_MW3.1_%02d.dat' % i)
            self.wave = de[:, 0] * 1e4
            dnue = np.abs(np.hstack([0, np.diff(3e18 / self.wave)]))
            for k, v in enumerate(np.arange(1, de.shape[1], 2)):
                norm = np.dot(dnue, de[:, v])
                DE1[j, k, :] = de[:, v] / norm
                norm = np.dot(dnue, de[:, v+1])
                DE2[j, k, :] = de[:, v+1] / norm
        shape = DE1.shape
        self.interpumin = LinearNDInterpolator(X, DE1.reshape(shape[0] *
                                                              shape[1],
                                                              shape[2]))
        self.interpumax = LinearNDInterpolator(X, DE2.reshape(shape[0] *
                                                              shape[1],
                                                              shape[2]))
    def get_nparams(self):
        ''' Return number of parameters '''
        if self.fixed:
            return 0
        else:
            return 3 

    def get_params(self):
        ''' Return current parameters '''
        if self.fixed:
            return []
        else:
            return [self.umin, self.gamma, self.qpah]

    def get_param_lims(self):
        ''' Return current parameter limits '''
        if self.fixed:
            return []
        else:
            return [self.umin_lims, self.gamma_lims, self.qpah_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        if self.fixed:
            return []
        else:
            return [self.umin_delta, self.gamma_delta, self.qpah_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        if self.fixed:
            return []
        else:
            return ['Umin', '$\gamma$', '$q_{pah}$']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        umin_flag = ((self.umin > self.umin_lims[0]) *
                     (self.umin < self.umin_lims[1]))
        gamma_flag = ((self.gamma > self.gamma_lims[0]) *
                      (self.gamma < self.gamma_lims[1]))
        qpah_flag = ((self.qpah > self.qpah_lims[0]) *
                     (self.qpah < self.qpah_lims[1]))
        return umin_flag * gamma_flag * qpah_flag

    def set_parameters_from_list(self, input_list, start_value):
        ''' Set parameters from a list and a start_value

        Parameters
        ----------
        input_list : list
            list of input parameters (could be much larger than number of
            parameters to be set)
        start_value : int
            initial index from list to read out parameters
        '''
        if not self.fixed:
            self.umin = input_list[start_value]
            self.gamma = input_list[start_value+1]
            self.qpah = input_list[start_value+2]

    def plot(self, ax, wave, color=[0/255., 175/255., 202/255.], alpha=0.2):
        ''' Plot Dust Law for given set of parameters '''
        dustem = self.evaluate(wave)
        ax.plot(wave, dustem, color=color, alpha=alpha)

    def evaluate(self, wave):
        ''' Evaluate Dust Law

        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength

        Returns
        -------
        DustE : numpy array (1 dim)
            Dust emission spectrum
        '''
        DustE = (self.interpumin(self.qpah, self.umin) * (1. - self.gamma) +
                 self.interpumax(self.qpah, self.umin) * self.gamma)

        return np.interp(wave, self.wave, DustE)
