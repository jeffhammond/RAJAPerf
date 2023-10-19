//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2017-23, Lawrence Livermore National Security, LLC
// and RAJA Performance Suite project contributors.
// See the RAJAPerf/LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

// Uncomment to add compiler directives for loop unrolling
//#define USE_RAJAPERF_UNROLL

#include "MASS3DEA.hpp"

#include "RAJA/RAJA.hpp"

#if defined(BUILD_STDPAR)

#include "common/StdParUtils.hpp"

#include <iostream>

namespace rajaperf {
namespace apps {

void MASS3DEA::runStdParVariant(VariantID vid,
                             size_t RAJAPERF_UNUSED_ARG(tune_idx)) {
#if defined(RUN_STDPAR)
  const Index_type run_reps = getRunReps();

  MASS3DEA_DATA_SETUP;

  switch (vid) {

  case Base_StdPar: {

    startTimer();
    for (RepIndex_type irep = 0; irep < run_reps; ++irep) {

      for (int e = 0; e < NE; ++e) {

        MASS3DEA_0_CPU

        CPU_FOREACH(d, x, MEA_D1D) {
          CPU_FOREACH(q, y, MEA_Q1D) {
            MASS3DEA_1
          }
        }

        MASS3DEA_2_CPU

        CPU_FOREACH(k1, x, MEA_Q1D) {
          CPU_FOREACH(k2, y, MEA_Q1D) {
            CPU_FOREACH(k3, z, MEA_Q1D) {
              MASS3DEA_3
            }
          }
        }

        CPU_FOREACH(i1, x, MEA_D1D) {
          CPU_FOREACH(i2, y, MEA_D1D) {
            CPU_FOREACH(i3, z, MEA_D1D) {
              MASS3DEA_4
            }
          }
        }

      } // element loop
    }
    stopTimer();

    break;
  }

  default:
    getCout() << "\n MASS3DEA : Unknown StdPar variant id = " << vid << std::endl;
  }

#else
  RAJA_UNUSED_VAR(vid);
#endif
}

} // end namespace apps
} // end namespace rajaperf

#endif  // BUILD_STDPAR

