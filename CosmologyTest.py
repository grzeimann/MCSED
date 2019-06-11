import numpy as np 
import matplotlib.pyplot as plt 
from cosmology import Cosmology as mcm
from astropy.cosmology import FlatLambdaCDM as acm
import timeit
from scipy.stats import linregress

tol = 1.0e-3
mc = mcm()
ac = acm(H0=69.,Om0=0.31)
zs = np.geomspace(0.0001,0.55,num=1001)
stepsize = np.geomspace(0.0000001,0.001,num=36)
best_step = stepsize[-1]*np.ones_like(zs) #Start out with largest stepsize possible (default)--will modify if necessary
acld = ac.luminosity_distance(zs).value
aclt = ac.lookback_time(zs).value

for i in range(len(zs)):
    if i%10==0: print "Starting routine for z =", zs[i]
    for j in range(len(stepsize)):
        mcld = mc.luminosity_distance(zs[i],stepsize=stepsize[j])/1.0e5
        err = (mcld-acld[i])/acld[i]
        if j==0: assert err<tol, "Error=%.3e"%(err)  #Make sure that we don't need to start out at even smaller stepsize
        if err >= tol: 
            best_step[i]=stepsize[j]
            if i%10==0: print best_step[i]
            break
a,b,r,p,stderr = linregress(np.log10(zs),np.log10(best_step))
print("slope: %f    intercept: %f" % (a, b))
print("R-squared: %f" % r**2)


plt.figure()
plt.loglog(zs,best_step)
plt.show()

# tic = timeit.default_timer()
# mcld = np.zeros_like(zs)
# mclt = np.zeros_like(zs)
# for i in range(len(zs)):
#     mcld[i] = mc.luminosity_distance(zs[i],stepsize=0.00001)
#     mclt[i] = mc.lookback_time(zs[i],stepsize=0.00001)
# toc = timeit.default_timer()
# mctime = toc-tic
# print "Timing for MCSED Luminosity Distance and Lookback Time:",toc-tic

# tic = timeit.default_timer()
# acld = ac.luminosity_distance(zs)
# aclt = ac.lookback_time(zs)
# toc = timeit.default_timer()
# actime = toc-tic
# print "Timing for Astropy Luminosity Distance and Lookback Time:",toc-tic

# plt.figure()
# plt.loglog(acld,mcld/1.0e5,'b.',label='')
# plt.loglog(acld,acld,'r--',label='1-1')
# plt.xlabel("Astropy Luminosity Distance (Mpc)")
# plt.ylabel("MCSED Luminosity Distance (Mpc)")
# plt.text(min(mcld)*1.1e-5,max(mcld)*0.1e-5,'MCSED Total Time: %.3f \nAstropy Total Time: %.3f'%(mctime,actime))
# plt.legend(loc='best')
# plt.savefig("LumDistCompStep00001.png")

# plt.figure()
# plt.loglog(aclt,mclt,'b.',label='')
# plt.loglog(aclt,aclt,'r--',label='1-1')
# plt.xlabel("Astropy Lookback Time (Gyr)")
# plt.ylabel("MCSED Lookback Time (Gyr)")
# plt.text(min(mclt)*1.1,max(mclt)*0.1,'MCSED Total Time: %.3f \nAstropy Total Time: %.3f'%(mctime,actime))
# plt.legend(loc='best')
# plt.savefig("LookTimeCompStep00001.png")