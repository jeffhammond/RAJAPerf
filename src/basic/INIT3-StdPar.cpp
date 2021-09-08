//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2017-21, Lawrence Livermore National Security, LLC
// and RAJA Performance Suite project contributors.
// See the RAJAPerf/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#include "INIT3.hpp"

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
namespace basic
{


void INIT3::runStdParVariant(VariantID vid)
{
#if defined(RUN_STDPAR)

  const Index_type run_reps = getRunReps();
  const Index_type ibegin = 0;
  const Index_type iend = getActualProblemSize();

#ifdef USE_RANGES
      auto range = std::views::iota(ibegin, iend);
      auto begin = std::begin(range);
      auto end   = std::end(range);
#else
      thrust::counting_iterator<Index_type> begin(ibegin);
      thrust::counting_iterator<Index_type> end(iend);
#endif

  INIT3_DATA_SETUP;

  auto init3_lam = [=](Index_type i) {
                     INIT3_BODY;
                   };

  switch ( vid ) {

    case Base_StdPar : {

      startTimer();
      for (RepIndex_type irep = 0; irep < run_reps; ++irep) {

        std::for_each( std::execution::par_unseq,
                       begin, end,
                       [=](Index_type i) {
          INIT3_BODY;
        });

      }
      stopTimer();

      break;
    }

    case Lambda_StdPar : {

      startTimer();
      for (RepIndex_type irep = 0; irep < run_reps; ++irep) {

        std::for_each( std::execution::par_unseq,
                       begin, end,
                       [=](Index_type i) {
          init3_lam(i);
        });

      }
      stopTimer();

      break;
    }

#if defined(RUN_RAJA_STDPAR)
    case RAJA_StdPar : {

      startTimer();
      for (RepIndex_type irep = 0; irep < run_reps; ++irep) {

        RAJA::forall<RAJA::simd_exec>(
          RAJA::RangeSegment(ibegin, iend), init3_lam);

      }
      stopTimer();

      break;
    }
#endif // RUN_RAJA_STDPAR

    default : {
      std::cout << "\n  INIT3 : Unknown variant id = " << vid << std::endl;
    }

  }

#endif
}

} // end namespace basic
} // end namespace rajaperf
