"""
Read in a series of FITS table files and make them accessible as
numarrays, optionally creating a HippoDraw NTuple.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
#
# $Id: FitsNTuple.py,v 1.8 2004/07/17 00:24:07 jchiang Exp $
#

import numpy as np
_cat = np.concatenate

class FitsNTupleError(RuntimeError):
    "Error reading in FITS table"

class FitsNTuple:
    def __init__(self, fitsfiles, extension=1):
        import sys, pyfits
        #
        # If fitsfile is not a list or tuple of file names, assume
        # it's a single file name and put it into a single element
        # tuple.
        #
        if type(fitsfiles) != type([]) and type(fitsfiles) != type(()):
            fitsfiles = (fitsfiles, )
        #
        # Process each file named in the list or tuple.
        #
        columnData = {}
        for i, file in enumerate(fitsfiles):
            #print "adding", file
            table = pyfits.open(file.strip(" "))
            if i == 0:
                self.names = table[extension].columns.names
            for name in self.names:
                myData = table[extension].data.field(name)
                if myData.dtype.name.find('float') == 0:
                    myData = np.array(myData, dtype=np.float)
                if myData.dtype.name.find('int') == 0:
                    myData = np.array(myData, dtype=np.int)
                if i == 0:
                    columnData[name] = myData
                else:
                    columnData[name] = _cat((columnData[name], myData))
        #
        # Add these columns to the internal dictionary.
        #
        self.__dict__.update(columnData)
    def makeNTuple(self, name=None, useNumArray=1):
        import hippo, sys, numarray
        if useNumArray:
            nt = hippo.NumArrayTuple()
        else:
            nt = hippo.NTuple()
        if name != None:
            nt.setTitle(name)
        ntc = hippo.NTupleController.instance()
        ntc.registerNTuple(nt)
        for name in self.names:
            if type(self.__dict__[name][0]) == np.ndarray:
                columns = self.__dict__[name]
                columns.transpose()
                for i, col in zip(xrange(sys.maxint), columns):
                    colname = "%s%i" % (name, i)
                    nt.addColumn(colname, col)
            else:
                nt.addColumn(name, self.__dict__[name])
        return nt
    def extend(self, other):
        for name in self.names:
            self.__dict__[name] = _cat((self.__dict__[name],
                                        other.__dict__[name]))
