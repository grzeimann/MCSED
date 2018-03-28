""" MCSED - sfh.py


1) Star Formation Histories:
    a) double power law:
    b) constant:

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
        self.nparams = 2

    def set_agelim(self, redshift):
        ''' Set the Age limit based on age of the universe '''
        C = Cosmology()
        self.age_lims[1] = np.log10(C.lookback_time(20.) -
                                    C.lookback_time(redshift))

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

    def plot(self, ax, color=[238/255., 90/255., 18/255.]):
        ''' Plot SFH for given set of parameters '''
        t = np.logspace(self.age_lims[0], self.age)
        sfr = self.evaluate(t)
        ax.plot(t, sfr, color=color, alpha=0.4)

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


class double_powerlaw:
    ''' The double powerlaw function provides a good description for the
    cosmic star formation history and any reasonable, smooth, continuous
    star formation rate.

        SFR(t) = 10**a * ((t/tau)**b + (t/tau)**-c)**-1
    '''
    def __init__(self, tau=-2.2, a=2., b=2., c=1., age=-.5, tau_lims=[-3., 1.],
                 a_lims=[-1., 5.], b_lims=[0., 5.], c_lims=[0., 5.],
                 age_lims=[-3., 0.4], tau_delta=0.2, a_delta=0.3, b_delta=0.3,
                 c_delta=0.3, age_delta=0.2):
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
        self.tau_delta = tau_delta
        self.a_delta = a_delta
        self.b_delta = b_delta
        self.c_delta = c_delta
        self.age_delta = age_delta
        self.age_lims = age_lims
        self.nparams = 5

    def set_agelim(self, redshift):
        ''' Set the Age limit based on age of the universe '''
        C = Cosmology()
        self.age_lims[1] = np.log10(C.lookback_time(20.) -
                                    C.lookback_time(redshift))

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

    def plot(self, ax, color=[238/255., 90/255., 18/255.]):
        ''' Plot SFH for given set of parameters '''
        t = np.logspace(self.age_lims[0], self.age)
        sfr = self.evaluate(t)
        ax.plot(t, sfr, color=color, alpha=0.4)

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
    def __init__(self, init_log_sfr=1.5, init_log_sfr_lims=[-5., 3.],
                 init_log_sfr_delta=0.2,
                 ages=[6.5, 7., 7.5, 8., 8.5, 9., 9.3]):
        ''' Initialize this class
        Parameters
        ----------
        TODO Fill these in
        '''
        self.ages = ages
        self.nparams = len(self.ages)
        self.nums = np.arange(1, self.nparams+1, dtype=int)
        for num in self.nums:
            setattr(self, 'sfr_' + str(num), init_log_sfr - num * 0.3)
            setattr(self, 'sfr_' + str(num) + '_lims', init_log_sfr_lims)
            setattr(self, 'sfr_' + str(num) + '_delta', init_log_sfr_delta)
        self.age_lims = [-3., self.ages[-1]-9.]

    def set_agelim(self, redshift):
        ''' Set the Age limit based on age of the universe '''
        C = Cosmology()
        self.age = np.log10(C.lookback_time(20.) -
                            C.lookback_time(redshift))

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

    def plot(self, ax, color=[238/255., 90/255., 18/255.]):
        ''' Plot SFH for given set of parameters '''
        t = np.array([6] + self.ages) - 9.
        sfr = self.evaluate(t)
        ax.step(10**t, np.hstack([sfr[0], sfr]), where='pre',
                color=color, alpha=0.4)

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
        return [10**p for p in self.get_params()]


class empirical:
    ''' The empirical SFH includes 6 bins of SFR at discrete time intervals '''
    def __init__(self, mass=9., age=-.5, mass_lims=[6., 12.], age_delta=0.3,
                 mass_delta=0.5, ages=[6.5, 7., 7.5, 8., 8.5, 9., 9.3]):
        ''' Initialize this class

        Parameters
        ----------
        TODO Fill these in
        '''
        self.ages = ages
        self.nparams = len(self.ages)+1
        self.nums = np.arange(1, self.nparams-1, dtype=int)
        self.mass = mass
        self.age = age
        self.age_delta = age_delta
        self.mass_lims = mass_lims
        self.mass_delta = mass_delta
        for num in np.arange(1,  self.nparams):
            setattr(self, 'frac_' + str(num), 1./(len(ages)))
            setattr(self, 'frac_' + str(num) + '_lims', [0., 1.])
            setattr(self, 'frac_' + str(num) + '_delta', 1./len(ages)/4.)
        self.age_lims = [self.ages[1]-9., self.ages[-1]-9.]

    def set_agelim(self, redshift):
        ''' Set the Age limit based on age of the universe '''
        C = Cosmology()
        self.age_lims[1] = np.log10(C.lookback_time(20.) -
                                    C.lookback_time(redshift))

    def get_params(self):
        ''' Return current parameters '''
        l = [self.mass, self.age]
        for num in self.nums:
            l.append(getattr(self, 'frac_' + str(num)))
        return l

    def get_param_lims(self):
        ''' Return current parameters limits '''
        l = [self.mass_lims, self.age_lims]
        for num in self.nums:
            l.append(getattr(self, 'frac_' + str(num) + '_lims'))
        return l

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        l = [self.mass_delta, self.age_delta]
        for num in self.nums:
            l.append(getattr(self, 'frac_' + str(num) + '_delta'))
        return l

    def get_names(self):
        ''' Return names of each parameter '''
        l = ['lmass', 'age']
        for num in self.nums:
            l.append('frac_' + str(num))
        return l

    def prior(self):
        ''' Uniform prior based on boundaries '''
        flag = True
        total = 0.
        for num in self.nums:
            val = getattr(self, 'frac_' + str(num))
            lims = getattr(self, 'frac_' + str(num) + '_lims')
            flag *= ((val > lims[0]) * (val < lims[1]))
            total += val
        mass_flag = ((self.mass > self.mass_lims[0]) *
                     (self.mass < self.mass_lims[1]))
        age_flag = ((self.age > self.age_lims[0]) *
                    (self.age < self.age_lims[1]))
        return flag * (total <= 1.) * mass_flag * age_flag

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
        self.age = input_list[start_value+1]
        total = 0
        for num in self.nums:
            v = input_list[start_value + num + 1]
            setattr(self, 'frac_' + str(num), v)
            total += v
        setattr(self, 'frac_' + str(self.nparams - 1), 1. - total)
        self.tdiff = np.diff(10**np.vstack([[0.] + self.ages,
                             [self.age+9.]*(len(self.ages)+1)]).min(axis=0))

    def plot(self, ax, color=[238/255., 90/255., 18/255.]):
        ''' Plot SFH for given set of parameters '''
        sel = np.where((np.array(self.ages)-9.) <= self.age)[0]
        t = np.hstack([6., np.array(self.ages)[sel]]) - 9
        sfr = np.array(self.evaluate(t))
        ax.step(10**t, np.hstack([sfr[0], sfr[sel]]), where='pre',
                color=color, alpha=0.2)

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
        v = self.get_params() + [getattr(self, 'frac_' +
                                               str(self.nparams - 1))]
        mass = 10**v[0]
        denominator = np.dot(self.tdiff, v[2:])

        return [mass * p / denominator for p in v[2:]]
