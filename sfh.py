""" MCSED - sfh.py


1) Star Formation Histories:
    a) double power law:
    b) constant:

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np


class double_powerlaw:
    ''' The double powerlaw function provides a good description for the
    cosmic star formation history and any reasonable, smooth, continuous
    star formation rate.

        SFR(t) = 10**a * ((t/tau)**b + (t/tau)**-c)**-1
    '''
    def __init__(self, tau=0., a=1., b=1., c=1., age=0., tau_lims=[-3., 1.],
                 a_lims=[-1., 3.], b_lims=[0., 5.], c_lims=[0., 5.],
                 age_lims=[0., 13.7]):
        ''' Initialize this class

        Input
        -----
        tau : float
            The inflection point between a rising SFR and a declining SFR
        a : float
            Normalization of the SFR(t) in log based 10
        b : float
            power for rising SFH
        c : float
            power for decling SFH
        tau_lims : list
            A two valued list for lower and upper boundary values for tau
        a_lims : list
            A two valued list for lower and upper boundary values for a
        b_lims : list
            A two valued list for lower and upper boundary values for b
        c_lims : list
            A two valued list for lower and upper boundary values for c
        '''
        self.tau = tau
        self.a = a
        self.b = b
        self.c = c
        self.age = age
        self.tau_lims = tau_lims
        self.a_lims = a_lims
        self.b_lims = b_lims
        self.c_lims = c_lims
        self.age_lims = age_lims
        self.nparams = 5

    def prior(self):
        ''' Uniform prior based on boundaries '''
        tau_flag = (self.tau > self.tau_lims[0])*(self.tau < self.tau_lims[1])
        a_flag = (self.a > self.a_lims[0])*(self.a < self.a_lims[1])
        b_flag = (self.b > self.b_lims[0])*(self.b < self.b_lims[1])
        c_flag = (self.c > self.c_lims[0])*(self.c < self.c_lims[1])
        age_flag = (self.age > self.age_lims[0])*(self.age < self.age_lims[1])
        return tau_flag * a_flag * b_flag * c_flag * age_flag

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
        self.a = input_list[start_value+1]
        self.b = input_list[start_value+2]
        self.c = input_list[start_value+3]
        self.age = input_list[start_value+4]

    def evaluate(self, t):
        ''' Evaluate double power law SFH

        Input
        -----
        t : numpy array (1 dim)
            time in Gyr

        Returns
        -------
        msfr : numpy array (1 dim)
            Star formation rate at given time in time array
        '''
        msfr = (10**(self.a) * ((t/self.tau)**self.b +
                                (t/self.tau)**(-self.c))**(-1))
        return msfr


class constant:
    def __init__(self, a=1.):
        self.a = a
        self.nparams = 1

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
        self.a = input_list[start_value]

    def evaluate(self, t):
        return self.a * np.ones(t.shape)
