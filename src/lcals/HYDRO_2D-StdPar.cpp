//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2017-21, Lawrence Livermore National Security, LLC
// and RAJA Performance Suite project contributors.
// See the RAJAPerf/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#include "HYDRO_2D.hpp"

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
namespace lcals
{


void HYDRO_2D::runStdParVariant(VariantID vid)
{
#if defined(RUN_STDPAR)

  const Index_type run_reps = getRunReps();
  const Index_type kbeg = 1;
  const Index_type kend = m_kn - 1;
  const Index_type jbeg = 1;
  const Index_type jend = m_jn - 1;

#ifdef USE_RANGES
  auto rangeK = std::views::iota(kbeg, kend);
  auto beginK = std::begin(rangeK);
  auto endK   = std::end(rangeK);
  auto rangeJ = std::views::iota(jbeg, jend);
  auto beginJ = std::begin(rangeJ);
  auto endJ   = std::end(rangeJ);
#else
  thrust::counting_iterator<Index_type> beginK(kbeg);
  thrust::counting_iterator<Index_type> endK(kend);
  thrust::counting_iterator<Index_type> beginJ(jbeg);
  thrust::counting_iterator<Index_type> endJ(jend);
#endif

  HYDRO_2D_DATA_SETUP;

  switch ( vid ) {

    case Base_StdPar : {

      startTimer();
      for (RepIndex_type irep = 0; irep < run_reps; ++irep) {

        std::for_each( std::execution::par,
                       beginK, endK,
                        [=](Index_type k) {
          std::for_each( std::execution::unseq,
                       beginJ, endJ,
                        [=](Index_type j) {
            HYDRO_2D_BODY1;
          });
        });

        std::for_each( std::execution::par,
                       beginK, endK,
                        [=](Index_type k) {
          std::for_each( std::execution::unseq,
                       beginJ, endJ,
                        [=](Index_type j) {
            HYDRO_2D_BODY2;
          });
        });

        std::for_each( std::execution::par,
                       beginK, endK,
                        [=](Index_type k) {
          std::for_each( std::execution::unseq,
                       beginJ, endJ,
                        [=](Index_type j) {
            HYDRO_2D_BODY3;
          });
        });

      }
      stopTimer();

      break;
    }

    case Lambda_StdPar : {

      auto hydro2d_base_lam1 = [=] (Index_type k, Index_type j) {
                                 HYDRO_2D_BODY1;
                               };
      auto hydro2d_base_lam2 = [=] (Index_type k, Index_type j) {
                                 HYDRO_2D_BODY2;
                               };
      auto hydro2d_base_lam3 = [=] (Index_type k, Index_type j) {
                                 HYDRO_2D_BODY3;
                               };

      startTimer();
      for (RepIndex_type irep = 0; irep < run_reps; ++irep) {

        std::for_each( std::execution::par,
                       beginK, endK,
                        [=](Index_type k) {
          std::for_each( std::execution::unseq,
                       beginJ, endJ,
                        [=](Index_type j) {
            hydro2d_base_lam1(k, j);
          });
        });

        std::for_each( std::execution::par,
                       beginK, endK,
                        [=](Index_type k) {
          std::for_each( std::execution::unseq,
                       beginJ, endJ,
                        [=](Index_type j) {
            hydro2d_base_lam2(k, j);
          });
        });

        std::for_each( std::execution::par,
                       beginK, endK,
                        [=](Index_type k) {
          std::for_each( std::execution::unseq,
                       beginJ, endJ,
                        [=](Index_type j) {
            hydro2d_base_lam3(k, j);
          });
        });

      }
      stopTimer();

      break;
    }

#if defined(RUN_RAJA_STDPAR)
    case RAJA_StdPar : {

      HYDRO_2D_VIEWS_RAJA;

      auto hydro2d_lam1 = [=] (Index_type k, Index_type j) {
                            HYDRO_2D_BODY1_RAJA;
                          };
      auto hydro2d_lam2 = [=] (Index_type k, Index_type j) {
                            HYDRO_2D_BODY2_RAJA;
                          };
      auto hydro2d_lam3 = [=] (Index_type k, Index_type j) {
                            HYDRO_2D_BODY3_RAJA;
                          };

      using EXECPOL =
        RAJA::KernelPolicy<
          RAJA::statement::For<0, RAJA::loop_exec,  // k
            RAJA::statement::For<1, RAJA::loop_exec,  // j
              RAJA::statement::Lambda<0>
            >
          >
        >;

      startTimer();
      for (RepIndex_type irep = 0; irep < run_reps; ++irep) {

        RAJA::kernel<EXECPOL>(
                     RAJA::make_tuple( RAJA::RangeSegment(kbeg, kend),
                                       RAJA::RangeSegment(jbeg, jend)),
                     hydro2d_lam1); 

        RAJA::kernel<EXECPOL>(
                     RAJA::make_tuple( RAJA::RangeSegment(kbeg, kend),
                                       RAJA::RangeSegment(jbeg, jend)),
                     hydro2d_lam2); 

        RAJA::kernel<EXECPOL>(
                     RAJA::make_tuple( RAJA::RangeSegment(kbeg, kend),
                                       RAJA::RangeSegment(jbeg, jend)),
                     hydro2d_lam3); 

      }
      stopTimer();

      break;
    }
#endif // RUN_RAJA_STDPAR

    default : {
      std::cout << "\n  HYDRO_2D : Unknown variant id = " << vid << std::endl;
    }

  }

#endif
}

} // end namespace lcals
} // end namespace rajaperf
