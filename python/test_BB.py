"""
@brief Comparison of C++ and pure Python implementations for the three
different modes (unbinned, binned, and point measurement).  The
reconstructions should be identical.

@author J. Chiang
"""
#
# $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/BayesianBlocks/python/test_BB.py,v 1.1.1.1 2011/09/03 00:55:59 jchiang Exp $
#
import numpy as num
import numpy.random as ra
import BayesianBlocks
import BayesianBlocks_python
import pylab_plotter as plot
from distributions import Gaussian
plot.pylab.ion()

BayesianBlocks.BayesianBlocks_enableFPE()

#
# Unbinned mode
#
interval0 = num.sort(ra.random(20))
interval1 = num.sort(ra.random(100)) + interval0[-1]
interval2 = num.sort(ra.random(20)) + interval1[-1]

seq = num.concatenate((interval0, interval1, interval2))
tmin = 0
tmax = seq[-1] + num.mean(seq[1:] - seq[:-1])

bb = BayesianBlocks.BayesianBlocks(seq, tmin, tmax)
bbp = BayesianBlocks_python.BayesianBlocks(seq)

xr = (0, 3)
bins = 50
binsize = float(xr[1] - xr[0])/bins
win0 = plot.Window(0)
plot.histogram(seq, xrange=xr, bins=bins)

for ncpPrior in range(1, 10):
    xx, yy = bb.lightCurve(ncpPrior)
    plot.curve(xx, num.array(yy)*binsize, color='r', linewidth=3)

    xxp, yyp = bbp.lightCurve(ncpPrior)
    plot.curve(xxp, num.array(yyp)*binsize, color='b', linewidth=1)

#
# Exposure weighting in unbinned mode
#
exposures = ra.random(len(seq))
bb.setCellSizes(exposures)

win1 = plot.Window(1)
for ncpPrior in range(1, 10):
    xx, yy = bb.lightCurve(ncpPrior)
    plot.curve(xx, num.array(yy)*binsize, color='r', linewidth=5)

#
# Binned mode
#
func0 = Gaussian(1, 10, 3)
class Histogram(object):
    def __init__(self, xmin, xmax, nx):
        self.bins = num.zeros(nx, dtype=float)
        self.xstep = (xmax - xmin)/float(nx-1)
        self.xvals = num.linspace(xmin, xmax, nx) + self.xstep/2.
    def add(self, x, wt=1):
        i = int((x - self.xvals[0])/self.xstep)
        if i < len(self.bins):
            self.bins[i] += wt

hist = Histogram(0, 20, 30)
for i in range(100):
    hist.add(func0.draw())

bb = BayesianBlocks.BayesianBlocks(hist.xvals[0] - hist.xstep/2.,
                                   hist.bins, 
                                   num.ones(len(hist.bins), dtype=float)*hist.xstep)
                                    
bbp = BayesianBlocks_python.BayesianBlocks(hist.bins, 
                                           num.ones(len(hist.bins))*hist.xstep,
                                           hist.xvals[0] - hist.xstep/2.)
win2 = plot.Window(2)
plot.xyplot(hist.xvals, hist.bins/hist.xstep, yerr=num.sqrt(hist.bins))
for ncpPrior in range(1, 10):
    xx, yy = bb.lightCurve(ncpPrior)
    plot.curve(xx, yy, color='r', linewidth=3)

    print (ncpPrior)
    xxp, yyp = bbp.lightCurve(ncpPrior)
    plot.curve(xxp, yyp, color='g', linewidth=1)

#
# Point measurement mode
#
func1 = Gaussian(1, 5, 1)
func2 = Gaussian(1, 10, 1)
func3 = Gaussian(1, 6, 1)

y = [func1.draw() for x in range(10)]
y.extend([func2.draw() for x in range(10)])
y.extend([func3.draw() for x in range(10)])
dy = num.ones(len(y))
x = range(len(y))

bb = BayesianBlocks.BayesianBlocks(x, y, dy)
bbp = BayesianBlocks_python.BayesianBlocks(x, y, dy)

win3 = plot.Window(3)

plot.xyplot(x, y, yerr=dy)
for ncpPrior in range(1, 10):
    xx, yy = bb.lightCurve(ncpPrior)
    plot.curve(xx, yy, color='r', linewidth=5)

    xxp, yyp = bbp.lightCurve(ncpPrior)
    plot.curve(xxp, yyp, color='g', linewidth=3)
