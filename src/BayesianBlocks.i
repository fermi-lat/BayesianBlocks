// -*- mode: c++ -*-
%module BayesianBlocks
%{
#include <fenv.h>
#include <utility>
#include <vector>
#include "BayesianBlocks/BayesianBlocks.h"
%}
%include stl.i
%exception {
   try {
      $action
   } catch (std::exception & eObj) {
      PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(eObj.what()));
      return NULL;
   }
}
%template(DoubleVector) std::vector<double>;
%template(PairDoubleVector) std::pair<std::vector<double>,std::vector<double> >;
%template(DoubleVectorVector) std::vector< std::vector<double> >;
%template(IntVector) std::vector<int>;
%include BayesianBlocks/BayesianBlocks.h
%extend BayesianBlocks {
   static void enableFPE() {
      feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
   }
   std::pair<std::vector<double>, std::vector<double> > 
      lightCurve(double ncp_prior) {
      std::vector<double> xx, yy;
      self->globalOpt(ncp_prior, xx, yy);
      return std::make_pair(xx, yy);
   }
}
