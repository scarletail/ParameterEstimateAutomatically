/* Simscape target specific file.
 * This file is generated for the Simscape network associated with the solver block 'RC1_Manual/Solver Configuration'.
 */

#include "ne_ds.h"
#include "RC1_Manual_c7e1c277_1_ds_obs_act.h"
#include "RC1_Manual_c7e1c277_1_ds_sys_struct.h"
#include "RC1_Manual_c7e1c277_1_ds_externals.h"
#include "RC1_Manual_c7e1c277_1_ds_external_struct.h"
#include "ssc_ml_fun.h"

int32_T RC1_Manual_c7e1c277_1_ds_obs_act(const NeDynamicSystem *sys, const
  NeDynamicSystemInput *t12, NeDsMethodOutput *t13)
{
  PmRealVector out;
  real_T nonscalar0[9];
  real_T nonscalar2[9];
  real_T Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C3;
  real_T Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C1;
  real_T zc_int1;
  boolean_T t3[2];
  real_T t5[1];
  real_T t6;
  real_T t8;
  ETTS0 efOut;
  size_t _in1ivar;
  real_T b_efOut[1];
  real_T U_idx_0;
  real_T X_idx_0;
  real_T X_idx_2;
  real_T X_idx_1;
  real_T X_idx_3;
  U_idx_0 = t12->mU.mX[0];
  X_idx_0 = t12->mX.mX[0];
  X_idx_1 = t12->mX.mX[1];
  X_idx_2 = t12->mX.mX[2];
  X_idx_3 = t12->mX.mX[3];
  out = t13->mOBS_ACT;
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
  Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C3 = -X_idx_0 +
    X_idx_2;
  Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C1 = X_idx_1 *
    -8.960573476702509E-6 + 1.0;
  t5[0ULL] = Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C1;
  t3[0ULL] = (Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C1 <
              0.1);
  t3[1ULL] = (Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C1 <=
              0.9);
  _in1ivar = 9ULL;
  tlu2_linear_nearest_prelookup((void *)&efOut.mField0, (void *)&efOut.mField1,
    (void *)&efOut.mField2, (void *)&efOut.mField3, (void *)nonscalar0, (void *)
    t5, (void *)t3, (void *)&_in1ivar);
  t3[0ULL] = (Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C1 <
              0.1);
  t3[1ULL] = (Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C1 <=
              0.9);
  _in1ivar = 9ULL;
  tlu2_1d_linear_nearest_value((void *)&b_efOut, (void *)efOut.mField0, (void *)
    efOut.mField1, (void *)efOut.mField2, (void *)efOut.mField3, (void *)
    nonscalar2, (void *)t3, (void *)&_in1ivar);
  zc_int1 = -(-U_idx_0 * b_efOut[0]) / -1.0;
  t6 = (zc_int1 -
        Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C3) / -1.0;
  t8 = -t6 / -1.0;
  out.mX[0] = t6;
  out.mX[1] = U_idx_0;
  out.mX[2] = U_idx_0;
  out.mX[3] = 0.0;
  out.mX[4] = -(-t6) / -1.0;
  out.mX[5] = 0.0;
  out.mX[6] = Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C1;
  out.mX[7] = -X_idx_3 + -U_idx_0;
  out.mX[8] = Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C3;
  out.mX[9] = X_idx_2;
  out.mX[10] = X_idx_0;
  out.mX[11] = X_idx_1 * 0.00027777777777777778;
  out.mX[12] = Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C1;
  out.mX[13] = U_idx_0;
  out.mX[14] = 0.0;
  out.mX[15] = X_idx_2;
  out.mX[16] = X_idx_2;
  out.mX[17] = Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C1;
  out.mX[18] = -U_idx_0;
  out.mX[19] = t6;
  out.mX[20] = Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C3;
  out.mX[21] = zc_int1;
  out.mX[22] = Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C1;
  out.mX[23] = X_idx_3;
  out.mX[24] = Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C3;
  out.mX[25] = X_idx_2;
  out.mX[26] = X_idx_0;
  out.mX[27] = Lithium_Cell_1RC_equivalent_circuit_model_single_temperature_C1;
  out.mX[28] = t6;
  out.mX[29] = 0.0;
  out.mX[30] = t8;
  out.mX[31] = 0.0;
  out.mX[32] = 0.0;
  out.mX[33] = t6;
  out.mX[34] = t8;
  (void)sys;
  (void)t13;
  return 0;
}
