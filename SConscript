# -*- python -*-
# @file SConscript
# @brief scons build specifications for BayesianBlocks
#
# $Header: /nfs/slac/g/glast/ground/cvs/ASP-scons/BayesianBlocks/SConscript,v 1.11 2011/08/12 19:15:21 jchiang Exp $
# Authors: J. Chiang <jchiang@slac.stanford.edu>
# Version: BayesianBlocks-04-00-00

import os
Import('baseEnv')
Import('listFiles')
libEnv = baseEnv.Clone()
swigEnv = baseEnv.Clone()

libEnv.Tool('BayesianBlocksLib', depsOnly = 1)
BayesianBlocksLib = libEnv.SharedLibrary('BayesianBlocks' ,
                                         ['src/BayesianBlocks.cxx'])

swigEnv.Tool('BayesianBlocksLib')
lib_BayesianBlocks = swigEnv.SwigLibrary('_BayesianBlocks', 
                                         'src/BayesianBlocks.i')

python_modules = listFiles(['python/*.py']) + ['src/BayesianBlocks.py']

swigEnv.Tool('registerTargets', package="BayesianBlocks",
             libraryCxts=[[BayesianBlocksLib, libEnv]],
             swigLibraryCxts=[[lib_BayesianBlocks, swigEnv]],
             includes=listFiles(['BayesianBlocks/*.h']),
             data = listFiles(['data/*']),
             python=python_modules)

