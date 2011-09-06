/**
 * @file BayesianBlocks
 * @brief Implementation of BB algorithm for event, binned and point
 * measurement data.
 * 
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/BayesianBlocks/BayesianBlocks/BayesianBlocks.h,v 1.1.1.1 2011/09/03 00:55:59 jchiang Exp $
 */

#ifndef _BayesianBlocks_h
#define _BayesianBlocks_h

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <deque>
#include <vector>

class BayesianBlocks {

public:

   BayesianBlocks(const std::vector<double> & arrival_times);

   BayesianBlocks(double start_time, 
                   const std::vector<double> & bin_content,
                   const std::vector<double> & bin_sizes);

   BayesianBlocks(const std::vector<double> & xx,
                   const std::vector<double> & yy,
                   const std::vector<double> & dy);

   /// @brief Compute the global optimum reconstruction of the 
   /// piecewise constant function.
   ///
   /// @param ncp_prior
   /// @param xx abscissa values of the reconstructed function
   /// @param yy ordinate values of the reconstructed function
   void globalOpt(double ncp_prior,
                  std::vector<double> & xx,
                  std::vector<double> & yy) const;

   typedef std::vector<double>::const_iterator const_iterator_t;

   double blockCost(size_t imin, size_t imax) const {
      return m_blockCost->operator()(imin, imax);
   }

   double blockSize(size_t imin, size_t imax) const;

   double blockContent(size_t imin, size_t imax) const;

   const std::vector<double> & cellContent() const {
      return m_cellContent;
   }

   const std::vector<double> & cellErrors() const {
      return m_cellErrors;
   }

   void setCellSizes(const std::vector<double> & cellSizes);

   /// @brief ncp_prior calibration for unbinned case as a function of
   /// number of events and false positive fraction.
   static double ncp_prior(double nevents, double fp_frac);

   static void enableFPE() {
#ifdef TRAP_FPE
      feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#else
      throw std::runtime_error("FPE handling is disabled "
                               "on non-linux platforms.");
#endif
   }

private:

   bool m_point_mode;
   bool m_binned;
   double m_tstart;
   std::vector<double> m_cellSizes;
   std::deque<double> m_unscaledCellSizePartialSums;
   std::deque<double> m_cellSizePartialSums;
   std::vector<double> m_cellContent;
   std::deque<double> m_cellContentPartialSums;
   std::vector<double> m_cellErrors;

   /// @brief Interface class for block cost functions.
   class BlockCost {
   public:

      BlockCost(const BayesianBlocks & bbObject) : m_bbObject(bbObject) {}

      virtual double operator()(size_t imin, size_t imax) const = 0;

   protected:

      const BayesianBlocks & m_bbObject;

   };

   /// @brief Cost function for unbinned or binned event-based data.
   class BlockCostEvent : public BlockCost {
   public:

      BlockCostEvent(const BayesianBlocks & bbObject) 
         : BlockCost(bbObject) {}

      virtual double operator()(size_t imin, size_t imax) const;

   };

   /// @brief Cost function for point measurement data.
   class BlockCostPoint : public BlockCost {
   public:

      BlockCostPoint(const BayesianBlocks & bbObject) 
         : BlockCost(bbObject) {}

      virtual double operator()(size_t imin, size_t imax) const;

   };

   BlockCost * m_blockCost;

   void generateCells(const std::vector<double> & arrival_times);
   void cellPartialSums();

   void ingestPointData(const std::vector<double> & xx,
                        const std::vector<double> & yy,
                        const std::vector<double> & dy);

   void lightCurve(const std::deque<size_t> & changePoints,
                   std::vector<double> & xx, 
                   std::vector<double> & yy) const;
};

#endif // _BayesianBlocks_h
