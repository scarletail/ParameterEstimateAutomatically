/* Simscape target specific file.
 * This file is generated for the Simscape network associated with the solver block 'RC1_Manual/Solver Configuration'.
 */

#include "ne_ds.h"
#include "RC1_Manual_c7e1c277_1_ds_dxf.h"
#include "RC1_Manual_c7e1c277_1_ds_sys_struct.h"
#include "RC1_Manual_c7e1c277_1_ds_externals.h"
#include "RC1_Manual_c7e1c277_1_ds_external_struct.h"
#include "ssc_ml_fun.h"

int32_T RC1_Manual_c7e1c277_1_ds_dxf(const NeDynamicSystem *sys, const
  NeDynamicSystemInput *t18, NeDsMethodOutput *t19)
{
  PmRealVector out;
  real_T nonscalar0[9];
  real_T nonscalar3[9];
  real_T nonscalar7[9];
  boolean_T t3[2];
  real_T t8[1];
  ETTS0 efOut;
  size_t _in1ivar;
  real_T b_efOut[1];
  real_T c_efOut[1];
  real_T d_efOut[1];
  real_T X_idx_1;
  real_T X_idx_3;
  X_idx_1 = t18->mX.mX[1];
  X_idx_3 = t18->mX.mX[3];
  out = t19->mDXF;
  nonscalar0[0] = 0.1;
  nonscalar0[1] = 0.2;
  nonscalar0[2] = 0.3;
  nonscalar0[3] = 0.4;
  nonscalar0[4] = 0.5;
  nonscalar0[5] = 0.6;
  nonscalar0[6] = 0.7;
  nonscalar0[7] = 0.8;
  nonscalar0[8] = 0.9;
  nonscalar3[0] = 0.000916949;
  nonscalar3[1] = 0.000627347;
  nonscalar3[2] = 0.000492104;
  nonscalar3[3] = 0.000421666;
  nonscalar3[4] = 0.000356116;
  nonscalar3[5] = 0.000291375;
  nonscalar3[6] = 0.00027227;
  nonscalar3[7] = 0.000262047;
  nonscalar3[8] = 0.000293451;
  nonscalar7[0] = 2.119;
  nonscalar7[1] = 2.143;
  nonscalar7[2] = 2.164;
  nonscalar7[3] = 2.185;
  nonscalar7[4] = 2.207;
  nonscalar7[5] = 2.235;
  nonscalar7[6] = 2.272;
  nonscalar7[7] = 2.324;
  nonscalar7[8] = 2.393;
  X_idx_1 = X_idx_1 * -8.960573476702509E-6 + 1.0;
  t8[0ULL] = X_idx_1;
  t3[0ULL] = (X_idx_1 < 0.1);
  t3[1ULL] = (X_idx_1 <= 0.9);
  _in1ivar = 9ULL;
  tlu2_linear_nearest_prelookup((void *)&efOut.mField0, (void *)&efOut.mField1,
    (void *)&efOut.mField2, (void *)&efOut.mField3, (void *)nonscalar0, (void *)
    t8, (void *)t3, (void *)&_in1ivar);
  t3[0ULL] = (X_idx_1 < 0.1);
  t3[1ULL] = (X_idx_1 <= 0.9);
  _in1ivar = 9ULL;
  tlu2_1d_linear_nearest_value((void *)&b_efOut, (void *)efOut.mField0, (void *)
    efOut.mField1, (void *)efOut.mField2, (void *)efOut.mField3, (void *)
    nonscalar3, (void *)t3, (void *)&_in1ivar);
  t3[0ULL] = (X_idx_1 < 0.1);
  t3[1ULL] = (X_idx_1 <= 0.9);
  _in1ivar = 9ULL;
  tlu2_1d_linear_nearest_derivatives((void *)&c_efOut, (void *)efOut.mField0,
    (void *)efOut.mField1, (void *)efOut.mField2, (void *)efOut.mField3, (void *)
    nonscalar3, (void *)t3, (void *)&_in1ivar);
  t3[0ULL] = (X_idx_1 < 0.1);
  t3[1ULL] = (X_idx_1 <= 0.9);
  _in1ivar = 9ULL;
  tlu2_1d_linear_nearest_derivatives((void *)&d_efOut, (void *)efOut.mField0,
    (void *)efOut.mField1, (void *)efOut.mField2, (void *)efOut.mField3, (void *)
    nonscalar7, (void *)t3, (void *)&_in1ivar);
  out.mX[0] = -(d_efOut[0] * -8.960573476702509E-6);
  out.mX[1] = -(X_idx_3 * (c_efOut[0] * -8.960573476702509E-6));
  out.mX[2] = -b_efOut[0];
  (void)sys;
  (void)t19;
  return 0;
}
