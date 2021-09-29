//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2017-21, Lawrence Livermore National Security, LLC
// and RAJA Performance Suite project contributors.
// See the RAJAPerf/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#include "POLYBENCH_3MM.hpp"

#include "RAJA/RAJA.hpp"

#ifdef USE_RANGES
#include <ranges>
#else
#include <thrust/iterator/counting_iterator.h>
#endif

#include <algorithm>
#include <execution>

#include <iostream>

#define USE_STDPAR_COLLAPSE 1

namespace rajaperf 
{
namespace polybench
{

void POLYBENCH_3MM::runStdParVariant(VariantID vid)
{
#if defined(RUN_STDPAR)

  const Index_type run_reps = getRunReps();

  POLYBENCH_3MM_DATA_SETUP;

#ifdef USE_RANGES
#  ifdef USE_STDPAR_COLLAPSE
  auto rangeIJ = std::views::iota((Index_type)0, ni*nj);
  auto beginIJ = std::begin(rangeIJ);
  auto endIJ   = std::end(rangeIJ);
  auto rangeIL = std::views::iota((Index_type)0, ni*nl);
  auto beginIL = std::begin(rangeIL);
  auto endIL   = std::end(rangeIL);
  auto rangeJL = std::views::iota((Index_type)0, nj*nl);
  auto beginJL = std::begin(rangeJL);
  auto endJL   = std::end(rangeJL);
#  else
  auto rangeI = std::views::iota((Index_type)0, ni);
  auto beginI = std::begin(rangeI);
  auto endI   = std::end(rangeI);
  auto rangeL = std::views::iota((Index_type)0, nl);
  auto beginL = std::begin(rangeL);
  auto endL   = std::end(rangeL);
#  endif
  auto rangeJ = std::views::iota((Index_type)0, nj);
  auto beginJ = std::begin(rangeJ);
  auto endJ   = std::end(rangeJ);
  auto rangeK = std::views::iota((Index_type)0, nk);
  auto beginK = std::begin(rangeK);
  auto endK   = std::end(rangeK);
  auto rangeM = std::views::iota((Index_type)0, nm);
  auto beginM = std::begin(rangeM);
  auto endM   = std::end(rangeM);
#else
#  ifdef USE_STDPAR_COLLAPSE
  thrust::counting_iterator<Index_type> beginIJ(0);
  thrust::counting_iterator<Index_type> endIJ(ni*nj);
  thrust::counting_iterator<Index_type> beginIL(0);
  thrust::counting_iterator<Index_type> endIL(ni*nl);
  thrust::counting_iterator<Index_type> beginJL(0);
  thrust::counting_iterator<Index_type> endJL(nj*nl);
#  else
  thrust::counting_iterator<Index_type> beginI(0);
  thrust::counting_iterator<Index_type> endI(ni);
  thrust::counting_iterator<Index_type> beginL(0);
  thrust::counting_iterator<Index_type> endL(nl);
#  endif
  thrust::counting_iterator<Index_type> beginJ(0);
  thrust::counting_iterator<Index_type> endJ(nj);
  thrust::counting_iterator<Index_type> beginK(0);
  thrust::counting_iterator<Index_type> endK(nk);
  thrust::counting_iterator<Index_type> beginM(0);
  thrust::counting_iterator<Index_type> endM(nm);
#endif

  switch ( vid ) {

    case Base_StdPar : {

      startTimer();
      for (RepIndex_type irep = 0; irep < run_reps; ++irep) {

#ifdef USE_STDPAR_COLLAPSE
        std::for_each( std::execution::par_unseq,
                       beginIJ, endIJ, [=](Index_type ij) {
            const auto i  = ij / nj;
            const auto j  = ij % nj;
#else
        std::for_each( std::execution::par_unseq,
                       beginI, endI, [=](Index_type i) {
          std::for_each(beginJ, endJ, [=](Index_type j) {
#endif
            POLYBENCH_3MM_BODY1;
            std::for_each(beginK, endK, [=,&dot](Index_type k) {
              POLYBENCH_3MM_BODY2;
            });
            POLYBENCH_3MM_BODY3;
#ifndef USE_STDPAR_COLLAPSE
          });
#endif
        });

#ifdef USE_STDPAR_COLLAPSE
        std::for_each( std::execution::par_unseq,
                       beginJL, endJL, [=](Index_type jl) {
            const auto j  = jl / nl;
            const auto l  = jl % nl;
#else
        std::for_each( std::execution::par_unseq,
                      beginJ, endJ, [=](Index_type j) {
          std::for_each(beginL, endL, [=](Index_type l) {
#endif
            POLYBENCH_3MM_BODY4;
            std::for_each(beginM, endM, [=,&dot](Index_type m) {
              POLYBENCH_3MM_BODY5;
            });
            POLYBENCH_3MM_BODY6;
#ifndef USE_STDPAR_COLLAPSE
          });
#endif
        });

#ifdef USE_STDPAR_COLLAPSE
        std::for_each( std::execution::par_unseq,
                       beginIL, endIL, [=](Index_type il) {
            const auto i  = il / nl;
            const auto l  = il % nl;
#else
        std::for_each( std::execution::par_unseq,
                       beginI, endI, [=](Index_type i) {
          std::for_each(beginL, endL, [=](Index_type l) {
#endif
            POLYBENCH_3MM_BODY7;
            std::for_each(beginJ, endJ, [=,&dot](Index_type j) {
              POLYBENCH_3MM_BODY8;
            });
            POLYBENCH_3MM_BODY9;
#ifndef USE_STDPAR_COLLAPSE
          });
#endif
        });

      }
      stopTimer();

      break;
    }

    case Lambda_StdPar : {

      auto poly_3mm_base_lam2 = [=] (Index_type i, Index_type j, Index_type k,
                                     Real_type &dot) {
                                  POLYBENCH_3MM_BODY2;
                                };
      auto poly_3mm_base_lam3 = [=] (Index_type i, Index_type j,
                                     Real_type &dot) {
                                  POLYBENCH_3MM_BODY3;
                                };
      auto poly_3mm_base_lam5 = [=] (Index_type j, Index_type l, Index_type m,
                                     Real_type &dot) {
                                   POLYBENCH_3MM_BODY5;
                                };
      auto poly_3mm_base_lam6 = [=] (Index_type j, Index_type l,
                                     Real_type &dot) {
                                  POLYBENCH_3MM_BODY6;
                                };
      auto poly_3mm_base_lam8 = [=] (Index_type i, Index_type l, Index_type j,
                                     Real_type &dot) {
                                  POLYBENCH_3MM_BODY8;
                                };
      auto poly_3mm_base_lam9 = [=] (Index_type i, Index_type l,
                                     Real_type &dot) {
                                  POLYBENCH_3MM_BODY9;
                                };

      startTimer();
      for (RepIndex_type irep = 0; irep < run_reps; ++irep) {

#ifdef USE_STDPAR_COLLAPSE
        std::for_each( std::execution::par_unseq,
                       beginIJ, endIJ, [=](Index_type ij) {
            const auto i  = ij / nj;
            const auto j  = ij % nj;
#else
        std::for_each( std::execution::par_unseq,
                       beginI, endI, [=](Index_type i) {
          std::for_each(beginJ, endJ, [=](Index_type j) {
#endif
            POLYBENCH_3MM_BODY1;
            std::for_each(beginK, endK, [=,&dot](Index_type k) {
              poly_3mm_base_lam2(i, j, k, dot);
            });
            poly_3mm_base_lam3(i, j, dot);
#ifndef USE_STDPAR_COLLAPSE
          });
#endif
        });

#ifdef USE_STDPAR_COLLAPSE
        std::for_each( std::execution::par_unseq,
                       beginJL, endJL, [=](Index_type jl) {
            const auto j  = jl / nl;
            const auto l  = jl % nl;
#else
        std::for_each( std::execution::par_unseq,
                       beginJ, endJ, [=](Index_type j) {
          std::for_each(beginL, endL, [=](Index_type l) {
#endif
            POLYBENCH_3MM_BODY4;
            std::for_each(beginM, endM, [=,&dot](Index_type m) {
              poly_3mm_base_lam5(j, l, m, dot);
            });
            poly_3mm_base_lam6(j, l, dot);
#ifndef USE_STDPAR_COLLAPSE
          });
#endif
        });

#ifdef USE_STDPAR_COLLAPSE
        std::for_each( std::execution::par_unseq,
                       beginIL, endIL, [=](Index_type il) {
            const auto i  = il / nl;
            const auto l  = il % nl;
#else
        std::for_each( std::execution::par_unseq,
                       beginI, endI, [=](Index_type i) {
          std::for_each(beginL, endL, [=](Index_type l) {
#endif
            POLYBENCH_3MM_BODY7;
            std::for_each(beginJ, endJ, [=,&dot](Index_type j) {
              poly_3mm_base_lam8(i, l, j, dot);
            });
            poly_3mm_base_lam9(i, l, dot);
#ifndef USE_STDPAR_COLLAPSE
          });
#endif
        });

      }
      stopTimer();

      break;
    }

#if defined(RUN_RAJA_STDPAR)
    case RAJA_StdPar : {

      POLYBENCH_3MM_VIEWS_RAJA;

      auto poly_3mm_lam1 = [=] (Real_type &dot) {
                                  POLYBENCH_3MM_BODY1_RAJA;
                                };
      auto poly_3mm_lam2 = [=] (Index_type i, Index_type j, Index_type k, 
                                Real_type &dot) {
                                  POLYBENCH_3MM_BODY2_RAJA;
                                };
      auto poly_3mm_lam3 = [=] (Index_type i, Index_type j,
                                Real_type &dot) {
                                  POLYBENCH_3MM_BODY3_RAJA;
                                };
      auto poly_3mm_lam4 = [=] (Real_type &dot) {
                                  POLYBENCH_3MM_BODY4_RAJA;
                                };
      auto poly_3mm_lam5 = [=] (Index_type j, Index_type l, Index_type m, 
                                Real_type &dot) {
                                  POLYBENCH_3MM_BODY5_RAJA;
                                };
      auto poly_3mm_lam6 = [=] (Index_type j, Index_type l,
                                Real_type &dot) {
                                  POLYBENCH_3MM_BODY6_RAJA;
                                };
      auto poly_3mm_lam7 = [=] (Real_type &dot) {
                                  POLYBENCH_3MM_BODY7_RAJA;
                                };
      auto poly_3mm_lam8 = [=] (Index_type i, Index_type l, Index_type j, 
                                Real_type &dot) {
                                  POLYBENCH_3MM_BODY8_RAJA;
                                };
      auto poly_3mm_lam9 = [=] (Index_type i, Index_type l,
                                Real_type &dot) {
                                  POLYBENCH_3MM_BODY9_RAJA;
                                };

      using EXEC_POL =
        RAJA::KernelPolicy<
          RAJA::statement::For<0, RAJA::loop_exec,
            RAJA::statement::For<1, RAJA::loop_exec,
              RAJA::statement::Lambda<0, RAJA::Params<0>>,
              RAJA::statement::For<2, RAJA::loop_exec,
                RAJA::statement::Lambda<1, RAJA::Segs<0,1,2>, RAJA::Params<0>>
              >,
              RAJA::statement::Lambda<2, RAJA::Segs<0,1>, RAJA::Params<0>>
            >
          >
        >;

      startTimer();
      for (RepIndex_type irep = 0; irep < run_reps; ++irep) {

        RAJA::kernel_param<EXEC_POL>(
          RAJA::make_tuple(RAJA::RangeSegment{0, ni},
                           RAJA::RangeSegment{0, nj},
                           RAJA::RangeSegment{0, nk}),
          RAJA::tuple<Real_type>{0.0},

          poly_3mm_lam1,
          poly_3mm_lam2,
          poly_3mm_lam3

        );

        RAJA::kernel_param<EXEC_POL>(
          RAJA::make_tuple(RAJA::RangeSegment{0, nj},
                           RAJA::RangeSegment{0, nl},
                           RAJA::RangeSegment{0, nm}),
          RAJA::tuple<Real_type>{0.0},

          poly_3mm_lam4,
          poly_3mm_lam5,
          poly_3mm_lam6

        ); 

        RAJA::kernel_param<EXEC_POL>(
          RAJA::make_tuple(RAJA::RangeSegment{0, ni},
                           RAJA::RangeSegment{0, nl},
                           RAJA::RangeSegment{0, nj}),
          RAJA::tuple<Real_type>{0.0},

          poly_3mm_lam7,
          poly_3mm_lam8,
          poly_3mm_lam9

        );

      } // end run_reps
      stopTimer();

      break;
    }
#endif // RUN_RAJA_STDPAR

    default : {
      std::cout << "\n  POLYBENCH_3MM : Unknown variant id = " << vid << std::endl;
    }

  }

#endif
}

} // end namespace polybench
} // end namespace rajaperf
