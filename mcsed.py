""" MCSED


1) CURRENT LIMITATIONS:
       A) Constant metallicity for input SSP
       B) Dust Emission is ad hoc from Draine and Li (2007)
   OPTIONAL FITTED PARAMETERS:
       A) SFH
           a) tau_sfh, age, a, b, c
       B) Dust law
           b) tau_dust, delta, Eb
   OUTPUT PRODUCTS:
       A) XXX Plot

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import logging

__all__ = ["Mcsed"]


class Mcsed:
    def __init__(self, filter_matrix, sfh_func=None, dust_abs_func=None,
                 data_mags=None, data_magerrs=None, redshift=None,
                 filter_flag=None):
        self.filter_matrix = filter_matrix
        self.sfh_func = sfh_func
        self.dust_abs_func = dust_abs_func
        self.data_mags = data_mags
        self.data_magerrs = data_magerrs
        self.redshift = redshift
        self.filter_flag = filter_flag

    def setup_logging(self):
        '''Setup Logging for MCSED

        Builds
        -------
        self.log : class
            self.log.info() is for general print and self.log.error() is
            for raise cases
        '''
        self.log = logging.getLogger('mcsed')
        if not len(self.log.handlers):
            # Set format for logger
            fmt = '[%(levelname)s - %(asctime)s] %(message)s'
            fmt = logging.Formatter(fmt)
            # Set level of logging
            level = logging.INFO
            # Set handler for logging
            handler = logging.StreamHandler()
            handler.setFormatter(fmt)
            handler.setLevel(level)
            # Build log with name, mcsed
            self.log = logging.getLogger('mcsed')
            # FIXME find out if this is needed
            self.log.setLevel(logging.DEBUG)
            self.log.addHandler(handler)
