"""
@brief Run BayesianBlocks analysis in unbinned mode on light curve
file processed by gtexposure.

@author J. Chiang
"""
#
# $Header$
#
import numpy as num
from BayesianBlocks import BayesianBlocks
from FitsNTuple import FitsNTuple

def bb_analysis(evfile, lcfile, fp_frac=1e-3, output='change_points.txt'):
    events = FitsNTuple(evfile)
    bb = BayesianBlocks(events.TIME)

    lc_nt = FitsNTuple(lcfile)
    bb.setCellSizes(lc_nt.EXPOSURE.tolist())

    nprior = BayesianBlocks.ncp_prior(len(events.TIME), fp_frac)
    x, y = bb.lightCurve(nprior)

    #
    # Compute fluxes over each block by summing exposures from light
    # curve file.
    #
    exposure = []
    counts = []
    xx = num.unique(x)

    for tmin, tmax in zip(xx[:-1], xx[1:]):
        imin = num.where(lc_nt.TIME > tmin)[0][0]
        imax = num.where(lc_nt.TIME < tmax)[0][-1] + 1
        exposure.append(sum(lc_nt.EXPOSURE[imin:imax]))
        counts.append(imax - imin)

    counts = num.array(counts)
    exposure = num.array(exposure)
    flux_vals = counts/exposure
    flux = []
    for item, yval in zip(flux_vals, y):
        flux.append(item)
        flux.append(item)

    output = open(output, 'w')
    for data in zip(x, y, flux):
        output.write('%.4f  %.4e  %.4e\n' % data)
    output.close()

    return len(num.unique(x)), nprior, len(events.TIME)

if __name__ == '__main__':
    import sys
    try:
        evfile, lcfile = sys.argv[1:3]
    except:
        print "usage: python bb_analysis.py evfile lcfile"
        sys.exit()
    ncp, ncpPrior, nevents = bb_analysis(fp_frac=fp_frac)
