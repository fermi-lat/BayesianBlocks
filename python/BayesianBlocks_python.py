"""
@brief Pure python implementation of the Bayesian Blocks algorithm
described by Jackson, Scargle et al. 2005, IEEE Signal Processing
Letters, 12, 105. (http://arxiv.org/abs/math/0309285)

@author J. Chiang <jchiang@slac.stanford.edu>
"""
#
# $Id: BayesianBlocks_python.py,v 1.1.1.1 2011/09/03 00:55:59 jchiang Exp $
#
import copy
import numpy as num

def gammln(xx):
    cof = [76.18009172947146, -86.50532032941677,
           24.01409824083091, -1.231739572450155,
           0.1208650973866179e-2, -0.5395239384953e-5]
    y = xx
    x = xx
    tmp = x + 5.5
    tmp -= (x + 0.5)*num.log(tmp)
    ser = 1.000000000190015
    for j in range(6):
        y += 1
        ser += cof[j]/y
    return -tmp + num.log(2.5066282746310005*ser/x)

class BayesianBlocks(object):
    """ 
    Unbinned mode:
    >>> bb = BayesianBlocks(arrival_times)

    Binned:
    >>> bb = BayesianBlocks(bin_content, bin_sizes, start_time)

    Point measurements:
    >>> bb = BayesianBlocks(time, flux, errors)

    Obtaining the piecewise constant light curve:
    >>> time, rate = bb.globalOpt(ncp_prior=1)
    """
    def __init__(self, *argv):
        self.point_mode = False
        self.use_ml = True
        if len(argv) == 1:
            events = list(argv[0])
            events.sort()
            events = num.array(events)
            self.cellContent = num.ones(len(argv[0]))
            self.cellSizes = self._generateCells(events)
            self.binned = False
        else:
            try:
                self._readPointData(argv)
            except TypeError:
                self.cellContent = copy.deepcopy(argv[0])
                self.cellSizes = copy.deepcopy(argv[1])
                self.tstart = argv[2]
            self.binned = True
    def _readPointData(self, argv):
        x, y, dy = (list(copy.deepcopy(argv[0])), 
                    list(copy.deepcopy(argv[1])),
                    list(copy.deepcopy(argv[2])))
        if len(x) != len(y) or len(y) != len(dy):
            raise RuntimeError("Point measurement mode: " + 
                               "input array sizes do not match")
        x.insert(0, x[0] - (x[1] - x[0]))
        x.append(x[-1] + (x[-1] - x[-2]))
        x = num.array(x)
        cell_bounds = (x[1:] + x[:-1])/2.
        self.tstart = cell_bounds[0]
        self.cellSizes = cell_bounds[1:] - cell_bounds[:-1]
        self.cellContent = y
        self.fluxes = num.array(y)
        self.errors = num.array(dy)
        self.point_mode = True
    def lightCurve(self, ncp_prior=1, use_ml=True):
        return self.globalOpt(ncp_prior, use_ml)
    def globalOpt(self, ncp_prior=1, use_ml=True):
        if self.point_mode:
            blockCost = self.blockCost_point
        else:
            blockCost = self.blockCost
        self.use_ml = use_ml
        opt, last = [], []
        opt.append(blockCost(0, 0) - ncp_prior)
        last.append(0)
        npts = len(self.cellContent)
        for nn in range(1, npts):
            max_opt = blockCost(0, nn) - ncp_prior
            jmax = 0
            for j in range(1, nn+1):
                my_opt = opt[j-1] + blockCost(j, nn) - ncp_prior
                if my_opt > max_opt:
                    max_opt = my_opt
                    jmax = j
            opt.append(max_opt)
            last.append(jmax)
        changePoints = []
        indx = last[-1]
        while indx > 0:
            changePoints.insert(0, indx)
            indx = last[indx-1]
        changePoints.insert(0, 0)
        changePoints.append(npts)
        return self._lightCurve(changePoints)
    def _lightCurve(self, changePoints):
        xx = []
        yy = []
        cell_sizes = self.cellSizes
        for imin, imax in zip(changePoints[:-1], changePoints[1:]):
            try:
                xx.extend([self.tstart + sum(cell_sizes[:imin]),
                           self.tstart + sum(cell_sizes[:imax])])
            except IndexError:
                xx.extend([self.tstart + imin*cell_sizes,
                           self.tstart + imax*cell_sizes])
            if self.point_mode:
                f, sig, weights = self._point_block_data(imin, imax-1)
                yval = sum(weights*f)
            else:
                yval = (sum(self.cellContent[imin:imax])
                        /sum(cell_sizes[imin:imax]))
            yy.extend([yval, yval])
        return xx, yy
    def _point_block_data(self, imin, imax):
        f, sig = self.fluxes[imin:imax+1], self.errors[imin:imax+1]
        weights = 1./sig**2/sum(1./sig**2)
        return f, sig, weights
    def blockCost_point(self, imin, imax):
        f, sig, weights = self._point_block_data(imin, imax)
        sigx2 = sum(weights*f**2) - (sum(weights*f))**2
        return -sigx2/2*sum(1./sig**2)
    def blockCost(self, imin, imax):
        size = self.blockSize(imin, imax)
        content = self.blockContent(imin, imax)
        if content == 0:
            return 0
        my_cost = content*(num.log(content/size) - 1)
        return my_cost
    def blockSize(self, imin, imax):
        try:
            return sum(self.cellSizes[imin:imax+1])
        except IndexError:
            return self.cellSizes*(imax - imin)
    def blockContent(self, imin, imax):
        return sum(self.cellContent[imin:imax+1])
    def _generateCells(self, events):
        self.tstart = (3*events[0] - events[1])/2.
        bounds = ((events[1:] + events[:-1])/2.).tolist()
        bounds.insert(0, self.tstart)
        bounds.append((3*events[-1] - events[-2])/2.)
        bounds = num.array(bounds)
        return bounds[1:] - bounds[:-1]

if __name__ == '__main__':
#    import hippoplotter as plot
#    import distributions as dist
#    nsamp = 200
#    events = dist.sample(dist.stepFunction(0.5, 0.7, amp=0.7), nsamp)
#
#    output = open('events.dat', 'w')
#    for event in events:
#        output.write("%12.4e\n" % event)
#    output.close()

    class Histogram(object):
        def __init__(self, xmin, xmax, nx):
            self.xmin = xmin
            self.dx = (xmax - xmin)/float(nx)
            self.binContent = num.zeros(nx)
            self.binSizes = self.dx*num.ones(nx)
        def add(self, xx, wt=1):
            indx = int((xx - self.xmin)/self.dx)
            self.binContent[indx] += wt

    events = [float(x.strip()) for x in open('events.dat', 'r')]

    hist = Histogram(0, 1, 50)
    for event in events:
        hist.add(event)

    bb = BayesianBlocks(events)
    xx, yy = bb.globalOpt(ncp_prior=1)

    bb2 = BayesianBlocks(hist.binContent, hist.binSizes, 0)
    xx2, yy2 = bb2.globalOpt(ncp_prior=1)

#    plot.histogram(events)
#    plot.scatter(xx, yy, oplot=1, pointRep='Line', color='red', autoscale=1)
#    plot.scatter(xx2, yy2, oplot=1, pointRep='Line', color='blue')
