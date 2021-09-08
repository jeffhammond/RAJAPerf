//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2017-21, Lawrence Livermore National Security, LLC
// and RAJA Performance Suite project contributors.
// See the RAJAPerf/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#include "FIR.hpp"

#include "RAJA/RAJA.hpp"

#ifdef USE_RANGES
#include <ranges>
#else
#include <thrust/iterator/counting_iterator.h>
#endif

#include <algorithm>
#include <execution>

#include <iostream>

namespace rajaperf 
{
namespace apps
{


void FIR::runStdParVariant(VariantID vid)
{
#if defined(RUN_STDPAR)

  const Index_type run_reps = getRunReps();
  const Index_type ibegin = 0;
  const Index_type iend = getActualProblemSize() - m_coefflen;

  FIR_COEFF;

  FIR_DATA_SETUP;

  Real_type coeff[FIR_COEFFLEN];
  std::copy(std::begin(coeff_array), std::end(coeff_array), std::begin(coeff));

  auto fir_lam = [=](Index_type i) {
                   FIR_BODY;
                 };
  
  switch ( vid ) {

    case Base_StdPar : {

#ifdef USE_RANGES
      auto range = std::views::iota(ibegin, iend);
      auto begin = std::begin(range);
      auto end   = std::end(range);
#else
      thrust::counting_iterator<Index_type> begin(ibegin);
      thrust::counting_iterator<Index_type> end(iend);
#endif

      startTimer();
      for (RepIndex_type irep = 0; irep < run_reps; ++irep) {

        std::for_each( std::execution::par_unseq,
                       begin, end,
                       [=](Index_type i) {
          FIR_BODY;
        });

      }
      stopTimer();

      break;
    } 

    case Lambda_StdPar : {

#ifdef USE_RANGES
      auto range = std::views::iota(ibegin,iend);
      auto begin = std::begin(range);
      auto end   = std::end(range);
#else
      thrust::counting_iterator<Index_type> begin(ibegin);
      thrust::counting_iterator<Index_type> end(iend);
#endif

      startTimer();
      for (RepIndex_type irep = 0; irep < run_reps; ++irep) {

        std::for_each( std::execution::par_unseq,
                       begin, end,
                       [=](Index_type i) {
           fir_lam(i);
        });

      }
      stopTimer();

      break;
    }

#if defined(RUN_RAJA_STDPAR)
    case RAJA_StdPar : {

      startTimer();
      for (RepIndex_type irep = 0; irep < run_reps; ++irep) {

        RAJA::forall<RAJA::loop_exec>(
          RAJA::RangeSegment(ibegin, iend), fir_lam);

      }
      stopTimer();

      break;
    }
#endif // RUN_RAJA_STDPAR

    default : {
      std::cout << "\n  FIR : Unknown variant id = " << vid << std::endl;
    }

  }

#endif
}

} // end namespace apps
} // end namespace rajaperf
