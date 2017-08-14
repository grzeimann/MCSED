""" MCSED - sfh.py


1) Star Formation Histories: 
    a) double power law:
    b) constant:

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np

def double_powerlaw(t, tau, a, b, c):
    msfr = 10**(a) * ( (t/tau)**b + (t/tau)**(-c))**(-1)
    return msfr 
    
def constant(t, a=1.):
    return a * np.ones(t.shape)