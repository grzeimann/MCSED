""" MCSED - sfh.py


1) Star Formation Histories:
    a) double power law:
    b) constant:

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np


class double_powerlaw:
    def __init__(self, tau=0., a=1., b=1., c=1.):
        self.tau = tau
        self.a = a
        self.b = b
        self.c = c
        self.nparams = 4

    def evaluate(self, t):
        msfr = (10**(self.a) * ((t/self.tau)**self.b +
                                (t/self.tau)**(-self.c))**(-1))
        return msfr


class constant:
    def __init__(self, a=1.):
        self.a = a
        self.nparams = 1

    def evaluate(self, t):
        return self.a * np.ones(t.shape)
