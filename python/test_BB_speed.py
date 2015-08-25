"""
@brief Plot execution times of C++ and pure Python implementations as a 
function of sample size.

@author J. Chiang
"""
#
# $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/BayesianBlocks/python/test_BB_speed.py,v 1.1.1.1 2011/09/03 00:55:59 jchiang Exp $
#
import sys
import time
import numpy as num
import numpy.random as ra
sys.path.insert(0, '/u/gl/jchiang/python')
from read_data import write_data

import BayesianBlocks
import BayesianBlocks_python

def running_time(BB_imp, nsamp=100, fp_frac=1e-3):
    events = ra.random(nsamp)
    events = num.sort(events)
    ncp_prior = BayesianBlocks.BayesianBlocks.ncp_prior(nsamp, fp_frac)
    tstart = time.time()
    dt_avg = num.mean(events[1:] - events[:-1])
    tmin = events[0] - dt_avg
    tmax = events[-1] + dt_avg
    bb = BB_imp.BayesianBlocks(events, tmin, tmax)
    x, y = bb.lightCurve(ncp_prior)
    dt = time.time() - tstart
    return dt, len(x)-2

if __name__ == '__main__':
    import pylab_plotter as plot
    plot.pylab.ion()

    nsamp = 200
    fp_frac = 1e-3
    nsamps = (10, 30, 100, 300, 1000)
    dts = [running_time(BayesianBlocks, nsamp=nsamp)[0] 
           for nsamp in nsamps]

#    nsamps = (10, 30, 100, 300, 1000)
#    dts_p = [running_time(BayesianBlocks_python, nsamp=nsamp)[0] 
#             for nsamp in nsamps[:4]]

    plot.xyplot(nsamps, dts, xlog=1, ylog=1)
#    plot.xyplot(nsamps[:4], dts_p, color='r',
#                yrange=(min(dts)/1.3, max(dts_p)*1.3),
#                xname='Sample size', yname='Running time (cpusec)')
