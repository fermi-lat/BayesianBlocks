"""
@brief Script to create light curve file, with exposure calculation,
used in BayesianBlocks analysis.

@author J. Chiang
"""
#
# $Header: /nfs/slac/g/glast/ground/cvs/BayesianBlocks/python/create_lc_file.py,v 1.3 2011/09/19 05:07:49 jchiang Exp $
#
import os
import numpy as num
import pyfits
from FitsNTuple import FitsNTuple
from GtApp import GtApp

gtbindef = GtApp('gtbindef')
gtbin = GtApp('gtbin')
gtexposure = GtApp('gtexposure')

def write_cell_boundaries(arrival_times, outfile, tbounds=None):
    times = num.array(arrival_times)
    cell_boundaries = ((times[:-1] + times[1:])/2.).tolist()
    if tbounds is None:
        tbounds = (3*times[0] - times[1])/2., (3*times[-1] - times[-2])/2.
    cell_boundaries.insert(0, tbounds[0])
    cell_boundaries.append(tbounds[1])
    output = open(outfile, 'w')
    for tmin, tmax in zip(cell_boundaries[:-1], cell_boundaries[1:]):
        output.write("%16.5f  %16.5f\n" % (tmin, tmax))
    output.close()

def tbounds(ft1file):
    ft1 = pyfits.open(ft1file)
    return ft1['EVENTS'].header['TSTART'], ft1['EVENTS'].header['TSTOP']

def ebounds(ft1file):
    ft1 = pyfits.open(ft1file)
    header = ft1['EVENTS'].header
    ndskeys = header['NDSKEYS']
    for i in range(1, ndskeys+1):
        type = 'DSTYP%i' % i
        if header[type] == 'ENERGY':
            key = 'DSVAL%i' % i
            data = header[key].split(':')
            return float(data[0]), float(data[1])

def create_lc_file(evfile, scfile, lcfile, irfs, specin=-2.,
                   tmp_ext='time_bins', clean=True, old_tbounds=False):
    binfile = '%s.txt' % tmp_ext
    tbinfile = '%s.fits' % tmp_ext
    events = FitsNTuple(evfile)
    tstart, tstop = tbounds(evfile)
    if old_tbounds:
        write_cell_boundaries(num.sort(events.TIME), binfile)
    else:
        write_cell_boundaries(num.sort(events.TIME), binfile,
                              tbounds=(tstart, tstop))
    gtbindef.run(bintype='t', binfile=binfile, outfile=tbinfile)
    emin, emax = ebounds(evfile)
    gtbin.run(evfile=evfile, scfile=scfile, outfile=lcfile, algorithm='LC',
              tbinalg='FILE', tstart=tstart, tstop=tstop, tbinfile=tbinfile)
    gtexposure.run(infile=lcfile, scfile=scfile, irfs=irfs, srcmdl='none',
                   specin=specin, emin=emin, emax=emax)
    if clean:
        os.remove(binfile)
        os.remove(tbinfile)

if __name__ == '__main__':
    import sys
    try:
        evfile, scfile, lcfile, irfs = sys.argv[1:5]
    except:
        print "usage: python create_lc_file.py evfile, scfile, lcfile, irfs"
        sys.exit()
    
    create_lc_file(evfile, scfile, lcfile, irfs)
