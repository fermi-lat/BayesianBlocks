/**
   @page BayesianBlocks
                  
   I implemented <a href="http://trotsky.arc.nasa.gov/~jeffrey/">Jeff
   Scargle</a>'s Bayesian Blocks algorithm for 1D data.  The main
   reference that I used is the article by <a
   href="http://trotsky.arc.nasa.gov/~jeffrey/paper_1d.pdf">Jackson,
   Scargle, et al. 2003</a>.  I also consulted Jeff's longer <a
   href="http://trotsky.arc.nasa.gov/~jeffrey/global.ps">Paper VI</a>
   for more details and the MatLab implementation.

   The basic method is intended to find the optimal partition of the
   interval such that the data in the interval can be described by a
   "piece-wise constant function", i.e., that the density of events in
   the interval can be modeled as a series of step functions.  For 1D
   Bayesian Blocks, each photon event is characterized by a single
   quantity, such as arrival time or measured energy.  For arrival
   times, the method gives an estimate of the light curve while for
   energies, it gives an estimate of the spectrum.  

*/
