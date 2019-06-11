""" MCSED - sfh.py


1) Star Formation Histories:
    a) constant
    b) exponential
    c) double powerlaw
    d) empirical

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np
from cosmology import Cosmology


class constant:
    ''' The constant star formation history '''
    def __init__(self, logsfr=1.0, age=-.5, logsfr_lims=[-3., 3.],
                 age_lims=[-3., 0.4], logsfr_delta=0.4, age_delta=0.2):
        ''' Initialize this class

        Parameters
        ----------
        logsfr : float
            Constant SFR in log based 10
        age : float
            Age of the galaxy when observed in log Gyrs
        logsfr_lims : list
            A two valued list for lower and upper boundary values for logsfr
        age_lims : list
            A two valued list for lower and upper boundary values for age
        logsfr_delta : float
            sigma to draw from a normal distribution when simulating galaxies
        age_delta : float
            sigma to draw from a normal distribution when simulating galaxies
        '''
        self.logsfr = logsfr
        self.age = age
        self.logsfr_lims = logsfr_lims
        self.age_lims = age_lims
        self.logsfr_delta = logsfr_delta
        self.age_delta = age_delta

    def set_agelim(self, redshift):
        ''' Set the Age limit based on age of the universe '''
        C = Cosmology()
        self.age_lims[1] = np.log10(C.lookback_time(20.) -
                                    C.lookback_time(redshift))

    def get_nparams(self):
        ''' Return number of parameters '''
        return 2

    def get_params(self):
        ''' Return current parameters '''
        return [self.logsfr, self.age]

    def get_param_lims(self):
        ''' Return current parameters '''
        return [self.logsfr_lims, self.age_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.logsfr_delta, self.age_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['Log SFR', 'Log Age']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        logsfr_flag = ((self.logsfr > self.logsfr_lims[0]) *
                       (self.logsfr < self.logsfr_lims[1]))
        age_flag = (self.age > self.age_lims[0])*(self.age < self.age_lims[1])
        return logsfr_flag * age_flag

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
        self.logsfr = input_list[start_value]
        self.age = input_list[start_value+1]

    def plot(self, ax, color=[238/255., 90/255., 18/255.], alpha=0.2):
        ''' Plot SFH for given set of parameters '''
        t = np.logspace(self.age_lims[0], self.age)
        sfr = self.evaluate(t)
        ax.plot(t, sfr, color=color, alpha=alpha)

    def evaluate(self, t):
        ''' Evaluate double power law SFH

        Parameters
        ----------
        t : numpy array (1 dim)
            time in Gyr

        Returns
        -------
        msfr : numpy array (1 dim)
            Star formation rate at given time in time array
        '''
        msfr = 10**self.logsfr * np.ones(t.shape)
        return msfr


class burst:
    ''' The constant star formation history '''
    def __init__(self, logsfr=1.0, age=-.5, burst_age=7.2, burst_sigma=0.4,
                 burst_strength=5., logsfr_lims=[-3., 3.], age_lims=[-3., 0.4],
                 burst_age_lims=[6., 7.5], burst_sigma_lims=[0.003, 0.5],
                 burst_strength_lims=[1., 10.], logsfr_delta=0.4,
                 age_delta=0.2, burst_age_delta=0.3, burst_sigma_delta=0.005,
                 burst_strength_delta=2.):
        ''' Initialize this class

        Parameters
        ----------
        logsfr : float
            Constant SFR in log based 10
        age : float
            Age of the galaxy when observed in log Gyrs
        logsfr_lims : list
        age_lims : list
            A two valued list for lower and upper boundary values for age
        logsfr_delta : float
            sigma to draw from a normal distribution when simulating galaxies
        age_delta : float
            sigma to draw from a normal distribution when simulating galaxies
        '''
        self.logsfr = logsfr
        self.age = age
        self.burst_age = burst_age
        self.burst_sigma = burst_sigma
        self.burst_strength = burst_strength
        self.logsfr_lims = logsfr_lims
        self.age_lims = age_lims
        self.burst_age_lims = burst_age_lims
        self.burst_sigma_lims = burst_sigma_lims
        self.burst_strength_lims = burst_strength_lims
        self.logsfr_delta = logsfr_delta
        self.age_delta = age_delta
        self.burst_age_delta = burst_age_delta
        self.burst_sigma_delta = burst_sigma_delta
        self.burst_strength_delta = burst_strength_delta

    def set_agelim(self, redshift):
        ''' Set the Age limit based on age of the universe '''
        C = Cosmology()
        self.age_lims[1] = np.log10(C.lookback_time(20.) -
                                    C.lookback_time(redshift))

    def get_nparams(self):
        ''' Return number of parameters '''
        return 4

    def get_params(self):
        ''' Return current parameters '''
        return [self.logsfr, self.age, self.burst_age,# self.burst_sigma,
                self.burst_strength]

    def get_param_lims(self):
        ''' Return current parameters '''
        return [self.logsfr_lims, self.age_lims, self.burst_age_lims,
                #self.burst_sigma_lims,
                self.burst_strength_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.logsfr_delta, self.age_delta, self.burst_age_delta,
                #self.burst_sigma_delta, 
                self.burst_strength_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['Log SFR', 'Log Age', 'Burst Age',# 'Burst Sigma',
                'Burst Strength']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        attr = ['logsfr', 'age', 'burst_age',# 'burst_sigma',
                'burst_strength']
        flag = True
        for att in attr:
            temp = getattr(self, att)
            templims = getattr(self, att + '_lims')
            flag *= (temp > templims[0]) * (temp < templims[1])
        return flag

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
        self.logsfr = input_list[start_value]
        self.age = input_list[start_value+1]
        self.burst_age = input_list[start_value+2]
        # self.burst_sigma = input_list[start_value+3]
        self.burst_strength = input_list[start_value+3]

    def plot(self, ax, color=[238/255., 90/255., 18/255.], alpha=0.2):
        ''' Plot SFH for given set of parameters '''
        t = np.logspace(self.age_lims[0], self.age)
        sfr = self.evaluate(t)
        ax.plot(t, sfr, color=color, alpha=alpha)

    def evaluate(self, t):
        ''' Evaluate double power law SFH

        Parameters
        ----------
        t : numpy array (1 dim)
            time in Gyr

        Returns
        -------
        msfr : numpy array (1 dim)
            Star formation rate at given time in time array
        '''
        nage = self.burst_age - 9.
        norm = (self.burst_strength * 10**self.logsfr /
                np.sqrt(2. * np.pi * self.burst_sigma))
        gauss = norm * np.exp(-0.5 * (np.log10(t) - nage)**2 / self.burst_sigma**2)
        msfr = 10**self.logsfr * np.ones(t.shape)
        return msfr + gauss


class polynomial:
    ''' The constant star formation history '''
    def __init__(self, age_locs=[6.5, 7.5, 8.5], age=-.5,
                 age_lims=[-3., 0.4], age_delta=0.2):
        ''' Initialize this class

        Parameters
        ----------
        logsfr : float
            Constant SFR in log based 10
        age : float
            Age of the galaxy when observed in log Gyrs
        logsfr_lims : list
            A two valued list for lower and upper boundary values for logsfr
        age_lims : list
            A two valued list for lower and upper boundary values for age
        logsfr_delta : float
            sigma to draw from a normal distribution when simulating galaxies
        age_delta : float
            sigma to draw from a normal distribution when simulating galaxies
        '''
        self.age_locs = age_locs
        self.nums = np.arange(1, len(age_locs)+1)
        for i in self.nums:
            setattr(self, 'p_' + str(i), 1.5 - i * 0.3)
            setattr(self, 'p_' + str(i) + '_lims', [-3., 3.])
            setattr(self, 'p_' + str(i) + '_delta', 0.3)
        self.middle_age = np.median(self.age_locs)
        x = np.array(self.age_locs) - self.middle_age
        L = []
        for i in self.nums:
            p = self.nums[-1] - i
            L.append(x**p)
        self.matrix = np.array(L).swapaxes(0, 1)
        self.age = age
        self.age_lims = age_lims
        self.age_delta = age_delta

    def set_agelim(self, redshift):
        ''' Set the Age limit based on age of the universe '''
        C = Cosmology()
        self.age_lims[1] = np.log10(C.lookback_time(20.) -
                                    C.lookback_time(redshift))
        self.age = self.age_lims[1]*1.

    def get_nparams(self):
        ''' Return number of parameters '''
        return len(self.age_locs) #+1

    def get_params(self):
        ''' Return current parameters '''
        p = []
        for i in self.nums:
            p.append(getattr(self, 'p_' + str(i)))
        #p.append(self.age)
        return p

    def get_param_lims(self):
        ''' Return current parameters '''
        p = []
        for i in self.nums:
            p.append(getattr(self, 'p_' + str(i) + '_lims'))
        #p.append(self.age_lims)
        return p

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        p = []
        for i in self.nums:
            p.append(getattr(self, 'p_' + str(i) + '_delta'))
        #p.append(self.age_delta)
        return p

    def get_names(self):
        ''' Return names of each parameter '''
        p = []
        for name in self.age_locs:
            name = 'P (%s)' % str(name)
            p.append(name)
        #p.append('Log Age')
        return p

    def prior(self):
        ''' Uniform prior based on boundaries '''
        flag = True
        for i in self.nums:
            p = getattr(self, 'p_' + str(i))
            p_lims = getattr(self, 'p_' + str(i) + '_lims')
            flag *= (p > p_lims[0]) * (p < p_lims[1])
        #age_flag = (self.age > self.age_lims[0])*(self.age < self.age_lims[1])
        return flag #* age_flag

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
        for i in self.nums:
            setattr(self, 'p_' + str(i), input_list[start_value + i - 1])
        #self.age = input_list[int(self.nums[-1]) + start_value]

    def plot(self, ax, color=[238/255., 90/255., 18/255.], alpha=0.2):
        ''' Plot SFH for given set of parameters '''
        t = np.logspace(self.age_lims[0], self.age)
        sfr = self.evaluate(t)
        ax.plot(t, sfr, color=color, alpha=alpha)

    def evaluate(self, t):
        ''' Evaluate double power law SFH

        Parameters
        ----------
        t : numpy array (1 dim)
            time in Gyr

        Returns
        -------
        msfr : numpy array (1 dim)
            Star formation rate at given time in time array
        '''
        v = np.array(self.get_params()[:])
        sol = np.linalg.lstsq(self.matrix, v)[0]
        msfr = 10**(np.polyval(sol, np.log10(t) - self.middle_age + 9.))
        return msfr


class exponential:
# WPBWPB comment correct? "constant" sfh?
    ''' The constant star formation history '''
    def __init__(self, logsfr=1.0, age=-1.0, tau=-1.5, logsfr_lims=[-3., 3.],
                 age_lims=[-3., 0.4], tau_lims=[-3.5, 3.5],
                 logsfr_delta=0.3, age_delta=0.2, tau_delta=0.4,
                 sign=1.0):
        ''' Initialize this class

        Parameters
        ----------
WPB REWRITE PARAMETERS
        logsfr : float
            Constant SFR in log based 10
        age : float
            Age of the galaxy when observed in log Gyrs
        logsfr_lims : list
            A two valued list for lower and upper boundary values for logsfr
        age_lims : list
            A two valued list for lower and upper boundary values for age
        logsfr_delta : float
            sigma to draw from a normal distribution when simulating galaxies
        age_delta : float
            sigma to draw from a normal distribution when simulating galaxies
        '''
        self.logsfr = logsfr
        self.age = age
        self.logsfr_lims = logsfr_lims
        self.age_lims = age_lims
        self.logsfr_delta = logsfr_delta
        self.age_delta = age_delta
        self.tau = tau
        self.tau_lims = tau_lims
        self.tau_delta = tau_delta
        self.sign = sign

    def set_agelim(self, redshift):
        ''' Set the Age limit based on age of the universe '''
        C = Cosmology()
        self.age_lims[1] = np.log10(C.lookback_time(20.) -
                                    C.lookback_time(redshift))

    def get_nparams(self):
        ''' Return number of parameters '''
        return 3

    def get_params(self):
        ''' Return current parameters '''
        return [self.logsfr, self.age, self.tau]

    def get_param_lims(self):
        ''' Return current parameters '''
        return [self.logsfr_lims, self.age_lims, self.tau_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.logsfr_delta, self.age_delta, self.tau_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['Log SFR', 'Log Age', r'Log $\tau$']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        logsfr_flag = ((self.logsfr > self.logsfr_lims[0]) *
                       (self.logsfr < self.logsfr_lims[1]))
        age_flag = (self.age > self.age_lims[0])*(self.age < self.age_lims[1])
        tau_flag = (self.tau > self.tau_lims[0])*(self.tau < self.age)
        return logsfr_flag * age_flag * tau_flag

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
        self.logsfr = input_list[start_value]
        self.age = input_list[start_value+1]
        self.tau = input_list[start_value+2]

    def plot(self, ax, color=[238/255., 90/255., 18/255.], alpha=0.2):
        ''' Plot SFH for given set of parameters '''
        t = np.logspace(self.age_lims[0], self.age)
        sfr = self.evaluate(t)
        ax.plot(t, sfr, color=color, alpha=alpha)

    def evaluate(self, t):
        ''' Evaluate double power law SFH

        Parameters
        ----------
        t : numpy array (1 dim)
            time in Gyr

        Returns
        -------
        msfr : numpy array (1 dim)
            Star formation rate at given time in time array
        '''
        if self.sign > 0.0:
            var = t
        else:
            var = 10**self.age - t
            var = np.max([var, np.zeros(var.shape)], axis=0)
        msfr = 10**self.logsfr * np.exp(-1. * var / 10**self.tau)
        return msfr


class double_powerlaw:
    ''' The double powerlaw function provides a good description for the
    cosmic star formation history and any reasonable, smooth, continuous
    star formation rate.

        SFR(t) = 10**a * ((t/tau)**b + (t/tau)**-c)**-1
    '''
    def __init__(self, tau=-2.4, a=3.0, b=2., c=1., age=-.5,
                 tau_lims=[-3., 1.], a_lims=[-1., 5.], b_lims=[0., 5.],
                 c_lims=[0., 5.], age_lims=[-3., 0.4], tau_delta=0.2,
                 a_delta=0.5, b_delta=0.5, c_delta=0.3, age_delta=0.2):
        ''' Initialize this class

        Parameters
        ----------
        tau : float
            The inflection point between a rising SFR and a declining SFR
        a : float
            Normalization of the SFR(t) in log based 10
        b : float
            power for rising SFH
        c : float
            power for decling SFH
        age : float
            Age of the galaxy when observed in log Gyrs
        age_lims : list
            A two valued list for lower and upper boundary values for age
        tau_lims : list
            A two valued list for lower and upper boundary values for tau
        a_lims : list
            A two valued list for lower and upper boundary values for a
        b_lims : list
            A two valued list for lower and upper boundary values for b
        c_lims : list
            A two valued list for lower and upper boundary values for c
        age_delta : float
            sigma to draw from a normal distribution when simulating galaxies
WPBWPB add: _delta for remaining parameters
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
        self.tau_delta = tau_delta
        self.a_delta = a_delta
        self.b_delta = b_delta
        self.c_delta = c_delta
        self.age_delta = age_delta
        self.age_lims = age_lims

    def set_agelim(self, redshift):
        ''' Set the Age limit based on age of the universe '''
        C = Cosmology()
        self.age_lims[1] = np.log10(C.lookback_time(20.) -
                                    C.lookback_time(redshift))

    def get_nparams(self):
        ''' Return number of parameters '''
        return 5

    def get_params(self):
        ''' Return current parameters '''
        return [self.tau, self.a, self.b, self.c, self.age]

    def get_param_lims(self):
        ''' Return current parameters '''
        return [self.tau_lims, self.a_lims, self.b_lims, self.c_lims,
                self.age_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.tau_delta, self.a_delta, self.b_delta, self.c_delta,
                self.age_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['$tau_{sfh}$', 'a', 'b', 'c', 'Log Age']

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

        Parameters
        ----------
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

    def plot(self, ax, color=[238/255., 90/255., 18/255.], alpha=0.2):
        ''' Plot SFH for given set of parameters '''
        t = np.logspace(self.age_lims[0], self.age)
        sfr = self.evaluate(t)
        ax.plot(t, sfr, color=color, alpha=alpha)

    def evaluate(self, t):
        ''' Evaluate double power law SFH

        Parameters
        ----------
        t : numpy array (1 dim)
            time in Gyr

        Returns
        -------
        msfr : numpy array (1 dim)
            Star formation rate at given time in time array
        '''
        t1 = 10**self.tau
        msfr = (10**(self.a) * ((t / t1)**self.b +
                                (t / t1)**(-self.c))**(-1))
        return msfr


class empirical_direct:
    ''' The empirical SFH includes 6 bins of SFR at discrete time intervals '''
    def __init__(self, init_log_sfr=1.2, init_log_sfr_lims=[-5., 3.],
                 init_log_sfr_delta=0.7,
                 ages=[8., 8.5, 9., 9.3]):
        ''' Initialize this class
        Parameters
        ----------
        TODO Fill these in
        '''
        self.ages = ages
        self.nums = np.arange(1, self.get_nparams()+1, dtype=int)
        for num in self.nums:
            setattr(self, 'sfr_' + str(num), init_log_sfr - num**1.2 * 0.3)
            setattr(self, 'sfr_' + str(num) + '_lims', init_log_sfr_lims)
            setattr(self, 'sfr_' + str(num) + '_delta', init_log_sfr_delta)
        self.sfr_1_delta = 0.3
        self.sfr_2_delta = 0.3
        self.age_lims = [-3., self.ages[-1]-9.]
        self.age = self.age_lims[1] * 1.

    def set_agelim(self, redshift):
        ''' Set the Age limit based on age of the universe '''
        C = Cosmology()
        self.age = np.log10(C.lookback_time(20.) -
                            C.lookback_time(redshift))

    def get_nparams(self):
        ''' Return number of parameters '''
        return len(self.ages)

    def get_params(self):
        ''' Return current parameters '''
        l = []
        for num in self.nums:
            l.append(getattr(self, 'sfr_' + str(num)))
        return l

    def get_param_lims(self):
        ''' Return current parameters limits '''
        l = []
        for num in self.nums:
            l.append(getattr(self, 'sfr_' + str(num) + '_lims'))
        return l

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        l = []
        for num in self.nums:
            l.append(getattr(self, 'sfr_' + str(num) + '_delta'))
        return l

    def get_names(self):
        ''' Return names of each parameter '''
        l = []
        for num in self.nums:
            l.append('sfr_' + str(num))
        return l

    def prior(self):
        ''' Uniform prior based on boundaries '''
        flag = True
        for num in self.nums:
            val = getattr(self, 'sfr_' + str(num))
            lims = getattr(self, 'sfr_' + str(num) + '_lims')
            flag *= ((val > lims[0]) * (val < lims[1]))
        return flag

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
        for num in self.nums:
            setattr(self, 'sfr_' + str(num), input_list[start_value + num - 1])

    def plot(self, ax, color=[238/255., 90/255., 18/255.], alpha=0.2):
        ''' Plot SFH for given set of parameters '''
        t = 10**(np.array([6] + self.ages) - 9.)
        sfr = self.evaluate(t)
#        ax.step(10**t, np.hstack([sfr[0], sfr]), where='pre',
#                color=color, alpha=alpha)
        ax.step(t, sfr, where='pre', color=color, alpha=alpha)


    def evaluate(self, t):
        ''' Evaluate empirical direct SFH
        Parameters
        ----------
        t : numpy array (1 dim)
            time in Gyr
        Returns
        -------
        msfr : numpy array (1 dim)
            Star formation rate at given time in time array
        '''
        if type(t) in [float, int]:
            t = np.array([t])
        elif type(t) != np.ndarray:
            t = np.array(t)
        # linear SFR in each SFH time bin
        sfr_bin = 10. ** np.array(self.get_params()) 
        # Ensure that self.ages, t are both in units of log years
        t_logyr = np.log10( t * 1e9 ) 
        bin_indx = np.searchsorted(self.ages, t_logyr, side="left")
        # adjust any times falling beyond the last SFH age bin
        bin_indx[ bin_indx > len(sfr_bin)-1 ] = len(sfr_bin)-1
        sfr = sfr_bin[ bin_indx ]
        return sfr


class empirical:
    ''' The empirical SFH includes 6 bins of SFR at discrete time intervals '''
    def __init__(self, mass=9.5, mass_lims=[6., 12.],
                 mass_delta=0.5, ages=[7., 8., 8.5, 9., 9.3]):
        ''' Initialize this class

        Parameters
        ----------
        TODO Fill these in
        '''
        self.ages = ages
        self.nums = np.arange(1, self.get_nparams(), dtype=int)
        self.mass = mass
        self.mass_lims = mass_lims
        self.mass_delta = mass_delta
        weight = 1. / (np.array(self.ages) - 6.5)**1.5
        weight /= weight.sum()
        for num, wei in zip(self.nums, weight[1:]):
            setattr(self, 'frac_' + str(num+1), wei)
            setattr(self, 'frac_' + str(num+1) + '_lims', [0., 1.])
            setattr(self, 'frac_' + str(num+1) + '_delta', wei/2.5)
        total = 0
        for num in self.nums:
            total += getattr(self, 'frac_' + str(num+1))
        self.firstbin = 1. - total
        self.age_lims = [self.ages[1]-9., self.ages[-1]-9.]
        self.age = self.age_lims[1] * 1.
        self.tdiff = np.diff(10**np.vstack([[0.] + self.ages,
                             [self.age+9.]*(len(self.ages)+1)]).min(axis=0))

    def set_agelim(self, redshift):
        ''' Set the Age limit based on age of the universe '''
        C = Cosmology()
        # self.age_lims[1] = np.log10(C.lookback_time(20.) -
        #                             C.lookback_time(redshift))
        self.age = np.log10(C.lookback_time(20.) - C.lookback_time(redshift))
        self.tdiff = np.diff(10**np.vstack([[0.] + self.ages,
                             [self.age+9.]*(len(self.ages)+1)]).min(axis=0))

    def get_nparams(self):
        ''' Return number of parameters '''
        return len(self.ages)

    def get_params(self):
        ''' Return current parameters '''
        l = [self.mass]
        for num in self.nums:
            l.append(getattr(self, 'frac_' + str(num+1)))
        return l

    def get_param_lims(self):
        ''' Return current parameters limits '''
        l = [self.mass_lims]
        for num in self.nums:
            l.append(getattr(self, 'frac_' + str(num+1) + '_lims'))
        return l

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        l = [self.mass_delta]
        for num in self.nums:
            l.append(getattr(self, 'frac_' + str(num+1) + '_delta'))
        return l

    def get_names(self):
        ''' Return names of each parameter '''
        l = ['lmass']
        for num in self.nums:
            l.append('frac_' + str(num+1))
        return l

    def prior(self):
        ''' Uniform prior based on boundaries '''
        flag = True
        total = 0.
        for num in self.nums:
            val = getattr(self, 'frac_' + str(num+1))
            lims = getattr(self, 'frac_' + str(num+1) + '_lims')
            flag *= ((val > lims[0]) * (val < lims[1]))
            total += val
        mass_flag = ((self.mass > self.mass_lims[0]) *
                     (self.mass < self.mass_lims[1]))
        return flag * (total <= 1.) * mass_flag

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
        self.mass = input_list[start_value]
        total = 0
        for num in self.nums:
            v = input_list[start_value + num]
            setattr(self, 'frac_' + str(num+1), v)
            total += v
        self.firstbin = 1. - total

    def plot(self, ax, color=[238/255., 90/255., 18/255.], alpha=0.2):
        ''' Plot SFH for given set of parameters '''
        sel = np.where((np.array(self.ages)-9.) <= self.age)[0]
        t = np.hstack([6., np.array(self.ages)[sel]]) - 9
        sfr = np.array(self.evaluate(t))
        ax.step(10**t, np.hstack([sfr[0], sfr[sel]]), where='pre',
                color=color, alpha=alpha)

    def evaluate(self, t):
        ''' Evaluate double power law SFH

        Parameters
        ----------
        t : numpy array (1 dim)
            time in Gyr

        Returns
        -------
        msfr : numpy array (1 dim)
            Star formation rate at given time in time array
        '''
        v = self.get_params()
        v.insert(1, self.firstbin)        
        mass = 10**v[0]
        denominator = np.dot(self.tdiff, v[1:])
        return [mass * p / denominator for p in v[1:]]
