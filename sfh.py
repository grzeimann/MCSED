""" MCSED - sfh.py


1) Star Formation Histories: 
    a) double power law:
    b) Noll:

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

def double_powerlaw(t, tau, a, b, c):
    msfr = 10**(a) * ( (t/tau)**b + (t/tau)**(-c))**(-1)
    return msfr 