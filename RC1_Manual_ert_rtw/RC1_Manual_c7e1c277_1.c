/* Simscape target specific file.
 * This file is generated for the Simscape network associated with the solver block 'RC1_Manual/Solver Configuration'.
 */

#include "ne_std.h"
#include "pm_default_allocator.h"
#include "ne_dae_fwd.h"
#include "ne_profiler_fwd.h"
#include "ne_dae_construct.h"
#include "nesl_la.h"
#include "RC1_Manual_c7e1c277_1_ds.h"

void RC1_Manual_c7e1c277_1_dae( NeDae **dae, const NeModelParameters
  *modelParams,
  const NeSolverParameters *solverParams)
{
  PmAllocator *ne_allocator;
  const McLinearAlgebra *linear_algebra_ptr =
    (solverParams->mLinearAlgebra == NE_FULL_LA) ?
    get_rtw_linear_algebra() :
    mc_get_csparse_linear_algebra();
  ne_allocator = pm_default_allocator();
  ne_dae_create( dae,
                RC1_Manual_c7e1c277_1_dae_ds( ne_allocator ),
                *solverParams,
                *modelParams,
                linear_algebra_ptr,
                NULL,
                NULL,
                NULL,
                ne_allocator);
}
