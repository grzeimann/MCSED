""" MCSED - dust.py


1) Dust Laws
    a) Calzetti:
    b) Noll:

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np


def calzettilaw(wave, Rv=4.05):
    ''' Calzetti et al. (2000) dust attenuation curve, k(wave)

    Parameters
    ----------
    wave : numpy array (1 dim)
        wavelength in Angstroms
    Rv : float
        extinction factor
    Returns
    -------
    k : numpy array (1 dim)
        A(wave) / Av = k(wave) / Rv
    '''
# WPB DELETE
#    print('this is Rv:  '+str(Rv))
    invwv = 1/(wave/1e4)
    sel1 = np.nonzero(wave < 0.63e4)[0]
    sel2 = np.nonzero(np.logical_and(wave >= 0.63e4, wave < 2.2e4))[0]
    k1 = np.zeros(sel1.shape)
    k2 = np.zeros(sel2.shape)
    k1 = (2.659 * (-2.156 + 1.509 * invwv[sel1] - 0.198 * invwv[sel1]**2 +
          0.011 * invwv[sel1]**3) + Rv)
    k2 = 2.659*(-1.857 + 1.040*invwv[sel2]) + Rv
    k = np.zeros(wave.shape)
    k[sel1] = k1
    k[sel2] = k2
    return k


class calzetti:
    ''' Calzetti Dust Law
    '''
    def __init__(self, EBV=0.15, EBV_lims=[-0.05, 1.50], EBV_delta=0.02, 
                 Rv=4.05, EBV_stars_gas=0.44):
        ''' Initialize Class

        Parameters
        -----
        EBV : float
            Color excess
            The observed (B-V) color minus the intrinsic (B-V) color
            Relates to Av, the magnitude of extinction in the V band (5500A) via
                Av = Rv * EBV
        Rv : float (held fixed throughout the fitting)
            Extinction factor
        EBV_stars_gas : float (held fixed throughout the fitting)
            coefficient between the attenuation applied to the stars and gas
            E(B-V)_stars = EBV_stars_gas * E(B-V)_gas
        '''
        self.EBV = EBV
        self.EBV_lims = EBV_lims
        self.EBV_delta = EBV_delta
        self.calz = None
        self.Rv = Rv
        self.EBV_stars_gas = EBV_stars_gas

    def get_nparams(self):
        ''' Return number of parameters '''
        return 1

    def get_params(self):
        ''' Return current parameters '''
        return [self.EBV]

    def get_param_lims(self):
        ''' Return current parameter limits '''
        return [self.EBV_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.EBV_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['E(B-V)']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        EBV_flag = (self.EBV > self.EBV_lims[0])*(self.EBV < self.EBV_lims[1])
        return EBV_flag

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
        self.EBV = input_list[start_value]

    def plot(self, ax, wave, color=[0/255., 175/255., 202/255.], alpha=0.2):
        ''' Plot Dust Law for given set of parameters '''
        dust = self.evaluate(wave)
        ax.plot(wave, dust, color=color, alpha=alpha)

    def evaluate(self, wave, new_wave=False):
        ''' Evaluate Dust Law

        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength
        new_wave : bool
            recompute k(wave), even if already present
            only used in emission line strengths: finer wavelength grid

        Returns
        -------
        Alam : numpy array (1 dim)
            Effective optical depth as a function of wavelength
            Observed = True * 10**(-0.4 * Av / Rv * k(wave))
            A(wave) = E(B-V) * k(wave) = Av / Rv * k(wave)
        '''
# WPB DELETE
#        print('this is self.Rv:  '+str(self.Rv))
        if self.calz is None:
            self.calz = calzettilaw(wave, self.Rv)
        if new_wave:
            kwave = calzettilaw(wave, self.Rv)
        else:
            kwave = self.calz
        Alam = self.EBV * kwave
        return Alam


class noll:
    ''' Prescription for dust law comes from Noll et al. (2009), with constants
    defined in Kriek & Conroy (2013).  This dust attenuation law includes a
    bump at 2175A and a modified Calzetti et al. (2000) attenuation curve.

        A(wave) = frac{A_V}{R_V} (k'(wave) + D(wave))
                     left(frac{wave}{5500}right) ^delta
        D(wave) = frac{E_b (wave,dellam)^2 }{(wave^2-lam0^2)^2
                     + (wave,dellam)^2}
    '''
    def __init__(self, EBV=0.15, delta=0.0, Eb=2.5, EBV_lims=[-0.05, 1.50],
                 delta_lims=[-1., 1.], Eb_lims=[-0.2, 6.], EBV_delta=0.02,
                 delta_delta=0.3, Eb_delta=1.0, 
                 Rv=4.05, EBV_stars_gas=0.44):
        ''' Initialize Class

        Parameters
        -----
        EBV : float
            Color excess
            The observed (B-V) color minus the intrinsic (B-V) color
            Relates to Av, the magnitude of extinction in the V band (5500A) via
                Av = Rv * EBV
        delta : float
            Power for powerlaw modification of Calzetti curve
        Eb : float
            Strength of 2175A bump.  See equation above for the Drude profile,
            D(wave)
        Rv : float (held fixed throughout the fitting)
            Extinction factor
        EBV_stars_gas : float (held fixed throughout the fitting)
            coefficient between the attenuation applied to the stars and gas
            E(B-V)_stars = EBV_stars_gas * E(B-V)_gas
        '''
        self.EBV = EBV
        self.delta = delta
        self.Eb = Eb
        self.EBV_lims = EBV_lims
        self.delta_lims = delta_lims
        self.Eb_lims = Eb_lims
        self.EBV_delta = EBV_delta
        self.delta_delta = delta_delta
        self.Eb_delta = Eb_delta
        self.calz = None
        self.Rv = Rv
        self.EBV_stars_gas = EBV_stars_gas

    def get_nparams(self):
        ''' Return number of parameters '''
        return 3

    def get_params(self):
        ''' Return current parameters '''
        return [self.EBV, self.delta, self.Eb]

    def get_param_lims(self):
        ''' Return current parameter limits '''
        return [self.EBV_lims, self.delta_lims, self.Eb_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.EBV_delta, self.delta_delta, self.Eb_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['E(B-V)', '$\delta$', '$E_b$']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        EBV_flag = (self.EBV > self.EBV_lims[0])*(self.EBV < self.EBV_lims[1])
        delta_flag = ((self.delta > self.delta_lims[0]) *
                      (self.delta < self.delta_lims[1]))
        Eb_flag = (self.Eb > self.Eb_lims[0])*(self.Eb < self.Eb_lims[1])
        return EBV_flag * delta_flag * Eb_flag

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
        self.EBV = input_list[start_value]
        self.delta = input_list[start_value+1]
        self.Eb = input_list[start_value+2]

    def plot(self, ax, wave, color=[0/255., 175/255., 202/255.], alpha=0.2):
        ''' Plot Dust Law for given set of parameters '''
        dust = self.evaluate(wave)
        ax.plot(wave, dust, color=color, alpha=alpha)

    def evaluate(self, wave, new_wave=False):
        ''' Evaluate Dust Law

        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength
        new_wave : bool
            recompute k(wave), even if already present
            only used in emission line strengths: finer wavelength grid

        Returns
        -------
        Alam : numpy array (1 dim)
            Effective optical depth as a function of wavelength
        '''
# WPB DELETE
#        print('this is self.Rv:  '+str(self.Rv))
        dellam = 350.
        lam0 = 2175.
        if self.calz is None:
            self.calz = calzettilaw(wave, self.Rv)
        if new_wave:
            kwave = calzettilaw(wave, self.Rv)
        else:
            kwave = self.calz

        Dlam = (self.Eb * (wave*dellam)**2 /
                          ((wave**2-lam0**2)**2+(wave*dellam)**2))
        Alam = (self.EBV * (kwave+Dlam)*(wave/5500)**(self.delta))
        return Alam


class noll_Eb_fixed:
    ''' Prescription for dust law comes from Noll et al. (2009), with constants
    defined in Kriek & Conroy (2013).  This dust attenuation law includes a
    bump at 2175A and a modified Calzetti et al. (2000) attenuation curve.

        A(wave) = frac{A_V}{R_V} (k'(wave) + D(wave))
                     left(frac{wave}{5500}right) ^delta
        D(wave) = frac{E_b (wave,dellam)^2 }{(wave^2-lam0^2)^2
                     + (wave,dellam)^2}

    In this variation, the dust bump is fixed at a particular value
    '''
    def __init__(self, EBV=0.15, delta=0.0, Eb=0., EBV_lims=[-0.05, 1.50],
                 delta_lims=[-1., 1.], EBV_delta=0.02, delta_delta=0.3,
                 Rv=4.05, EBV_stars_gas=0.44):
        ''' Initialize Class

        Parameters
        -----
        EBV : float
            Color excess
            The observed (B-V) color minus the intrinsic (B-V) color
            Relates to Av, the magnitude of extinction in the V band (5500A) via
                Av = Rv * EBV
        delta : float
            Power for powerlaw modification of Calzetti curve
        Eb : float (held fixed throughout the fitting)
            Strength of 2175A bump.  See equation above for the Drude profile,
            D(wave)
        Rv : float (held fixed throughout the fitting)
            Extinction factor
        EBV_stars_gas : float (held fixed throughout the fitting)
            coefficient between the attenuation applied to the stars and gas
            E(B-V)_stars = EBV_stars_gas * E(B-V)_gas
        '''
        self.EBV = EBV
        self.delta = delta
        self.Eb = Eb
        self.EBV_lims = EBV_lims
        self.delta_lims = delta_lims
        self.EBV_delta = EBV_delta
        self.delta_delta = delta_delta
        self.calz = None
        self.Rv = Rv
        self.EBV_stars_gas = EBV_stars_gas

    def get_nparams(self):
        ''' Return number of parameters '''
        return 2

    def get_params(self):
        ''' Return current parameters '''
        return [self.EBV, self.delta]

    def get_param_lims(self):
        ''' Return current parameter limits '''
        return [self.EBV_lims, self.delta_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.EBV_delta, self.delta_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['E(B-V)', '$\delta$']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        EBV_flag = (self.EBV > self.EBV_lims[0])*(self.EBV < self.EBV_lims[1])
        delta_flag = ((self.delta > self.delta_lims[0]) *
                      (self.delta < self.delta_lims[1]))
        return EBV_flag * delta_flag

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
        self.EBV = input_list[start_value]
        self.delta = input_list[start_value+1]

    def plot(self, ax, wave, color=[0/255., 175/255., 202/255.], alpha=0.2):
        ''' Plot Dust Law for given set of parameters '''
        dust = self.evaluate(wave)
        ax.plot(wave, dust, color=color, alpha=alpha)

    def evaluate(self, wave, new_wave=False):
        ''' Evaluate Dust Law

        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength
        new_wave : bool
            recompute k(wave), even if already present
            only used in emission line strengths: finer wavelength grid

        Returns
        -------
        Alam : numpy array (1 dim)
            Effective optical depth as a function of wavelength
        '''
        dellam = 350.
        lam0 = 2175.
        if self.calz is None:
            self.calz = calzettilaw(wave, self.Rv)
        if new_wave:
            kwave = calzettilaw(wave, self.Rv)
        else:
            kwave = self.calz

        Dlam = (self.Eb * (wave*dellam)**2 /
                          ((wave**2-lam0**2)**2+(wave*dellam)**2))
        Alam = (self.EBV * (kwave+Dlam)*(wave/5500)**(self.delta))
        return Alam



class reddy:
    ''' 
    WPBWPB fill from mallory's description
    '''
    def __init__(self, EBV=0.15, EBV_lims=[-0.15, 2.40], EBV_delta=0.02, 
                 Rv=2.505, EBV_stars_gas=0.44):
        ''' Initialize Class

        Parameters
        -----
        EBV : float
            Color excess
            The observed (B-V) color minus the intrinsic (B-V) color
            Relates to Av, the magnitude of extinction in the V band (5500A) via
                Av = Rv * EBV
        Rv : float (held fixed throughout the fitting)
            Extinction factor
        EBV_stars_gas : float (held fixed throughout the fitting)
            coefficient between the attenuation applied to the stars and gas
            E(B-V)_stars = EBV_stars_gas * E(B-V)_gas
        '''
        self.EBV = EBV
        self.EBV_lims = EBV_lims
        self.EBV_delta = EBV_delta
        self.klam = None
        self.Rv = Rv
        self.EBV_stars_gas = EBV_stars_gas

    def get_nparams(self):
        ''' Return number of parameters '''
        return 1

    def get_params(self):
        ''' Return current parameters '''
        return [self.EBV]

    def get_param_lims(self):
        ''' Return current parameter limits '''
        return [self.EBV_lims]

    def get_param_deltas(self):
        ''' Return current parameter deltas '''
        return [self.EBV_delta]

    def get_names(self):
        ''' Return names of each parameter '''
        return ['E(B-V)']

    def prior(self):
        ''' Uniform prior based on boundaries '''
        EBV_flag = (self.EBV > self.EBV_lims[0])*(self.EBV < self.EBV_lims[1])
        return EBV_flag

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
        self.EBV = input_list[start_value]

    def plot(self, ax, wave, color=[0/255., 175/255., 202/255.], alpha=0.2):
        ''' Plot Dust Law for given set of parameters '''
        dust = self.evaluate(wave)
        ax.plot(wave, dust, color=color, alpha=alpha)

    def reddylaw(self, wave):
        ''' Reddy et al. (2015) dust attenuation curve, k(wave)
    
        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength in Angstroms
        Returns
        -------
        k : numpy array (1 dim)
            A(wave) / Av = k(wave) / Rv
        '''
# WPB DELETE
#        print('this is self.Rv:  '+str(self.Rv))
        Rv = self.Rv
        invwv = 1/(wave/1e4)
        sel1 = np.nonzero(wave < 0.60e4)[0]
        sel2 = np.nonzero(np.logical_and(wave >= 0.60e4, wave < 2.2e4))[0]
        k1 = np.zeros(sel1.shape)
        k2 = np.zeros(sel2.shape)
        k1 = (-5.726 + 4.004 * invwv[sel1] - 0.525 * invwv[sel1]**2 +
              0.029 * invwv[sel1]**3 + Rv)
        k2 = (-2.672 - 0.010 * invwv[sel2] + 1.532 * invwv[sel1]**2 -
              0.412 * invwv[sel1]**3 + Rv)
        k = np.zeros(wave.shape)
        k[sel1] = k1
        k[sel2] = k2
        return k

    def evaluate(self, wave, new_wave=False):
        ''' Evaluate Dust Law

        Parameters
        ----------
        wave : numpy array (1 dim)
            wavelength
        new_wave : bool
            recompute k(wave), even if already present
            only used in emission line strengths: finer wavelength grid

        Returns
        -------
        Alam : numpy array (1 dim)
            Effective optical depth as a function of wavelength
            Observed = True * 10**(-0.4 * Av / Rv * k(wave))
            A(wave) = E(B-V) * k(wave) = Av / Rv * k(wave)
        '''
        if self.klam is None:
            #WPBWPB delete
            print('self.klam not yet defined')
            self.klam = self.reddylaw(wave)
        else:
            print('self.klam already defined')
        if new_wave:
            kwave = reddylaw(wave, self.Rv)
        else:
            kwave = self.klam
        Alam = self.EBV * kwave
        return Alam

