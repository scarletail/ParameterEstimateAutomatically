/* Simscape target specific file.
 * This file is generated for the Simscape network associated with the solver block 'RC1_Manual/Solver Configuration'.
 */

#include "ne_ds.h"
#include "RC1_Manual_c7e1c277_1_ds_duy.h"
#include "RC1_Manual_c7e1c277_1_ds_sys_struct.h"
#include "RC1_Manual_c7e1c277_1_ds_externals.h"
#include "RC1_Manual_c7e1c277_1_ds_external_struct.h"
#include "ssc_ml_fun.h"

int32_T RC1_Manual_c7e1c277_1_ds_duy(const NeDynamicSystem *sys, const
  NeDynamicSystemInput *t8, NeDsMethodOutput *t9)
{
  PmRealVector out;
  real_T nonscalar0[9];
  real_T nonscalar2[9];
  boolean_T t3[2];
  real_T t4[1];
  ETTS0 efOut;
  size_t _in1ivar;
  real_T b_efOut[1];
  real_T X_idx_1;
  X_idx_1 = t8->mX.mX[1];
  out = t9->mDUY;
  nonscalar0[0] = 0.1;
  nonscalar0[1] = 0.2;
  nonscalar0[2] = 0.3;
  nonscalar0[3] = 0.4;
  nonscalar0[4] = 0.5;
  nonscalar0[5] = 0.6;
  nonscalar0[6] = 0.7;
  nonscalar0[7] = 0.8;
  nonscalar0[8] = 0.9;
  nonscalar2[0] = 0.00064444;
  nonscalar2[1] = 0.000570341;
  nonscalar2[2] = 0.000518561;
  nonscalar2[3] = 0.000503711;
  nonscalar2[4] = 0.000474134;
  nonscalar2[5] = 0.000474071;
  nonscalar2[6] = 0.000459256;
  nonscalar2[7] = 0.000459256;
  nonscalar2[8] = 0.000459266;
  X_idx_1 = X_idx_1 * -8.960573476702509E-6 + 1.0;
  t4[0ULL] = X_idx_1;
  t3[0ULL] = (X_idx_1 < 0.1);
  t3[1ULL] = (X_idx_1 <= 0.9);
  _in1ivar = 9ULL;
  tlu2_linear_nearest_prelookup((void *)&efOut.mField0, (void *)&efOut.mField1,
    (void *)&efOut.mField2, (void *)&efOut.mField3, (void *)nonscalar0, (void *)
    t4, (void *)t3, (void *)&_in1ivar);
  t3[0ULL] = (X_idx_1 < 0.1);
  t3[1ULL] = (X_idx_1 <= 0.9);
  _in1ivar = 9ULL;
  tlu2_1d_linear_nearest_value((void *)&b_efOut, (void *)efOut.mField0, (void *)
    efOut.mField1, (void *)efOut.mField2, (void *)efOut.mField3, (void *)
    nonscalar2, (void *)t3, (void *)&_in1ivar);
  out.mX[0] = -(-(-(-(-(-b_efOut[0])))));
  (void)sys;
  (void)t9;
  return 0;
}
