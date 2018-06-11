/*
 * File: RC1_Manual.h
 *
 * Code generated for Simulink model 'RC1_Manual'.
 *
 * Model version                  : 1.32
 * Simulink Coder version         : 8.14 (R2018a) 06-Feb-2018
 * C/C++ source code generated on : Mon Jun 11 10:34:03 2018
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_RC1_Manual_h_
#define RTW_HEADER_RC1_Manual_h_
#include <string.h>
#include <math.h>
#ifndef RC1_Manual_COMMON_INCLUDES_
# define RC1_Manual_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "nesl_rtw.h"
#include "RC1_Manual_c7e1c277_1_gateway.h"
#endif                                 /* RC1_Manual_COMMON_INCLUDES_ */

#include "math.h"
#include "rt_matrixlib.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

/* Forward declaration for rtModel */
typedef struct tag_RTM RT_MODEL;

/* Block signals and states (default storage) for system '<Root>' */
typedef struct {
  real_T STATE_1[5];                   /* '<S9>/STATE_1' */
  real_T INPUT_1_1_1[4];               /* '<S9>/INPUT_1_1_1' */
  real_T INPUT_1_1_1_discrete[2];      /* '<S9>/INPUT_1_1_1' */
  real_T OUTPUT_1_1;                   /* '<S9>/OUTPUT_1_1' */
  void* STATE_1_Simulator;             /* '<S9>/STATE_1' */
  void* STATE_1_SimulationData;        /* '<S9>/STATE_1' */
  void* STATE_1_DiagnosticManager;     /* '<S9>/STATE_1' */
  void* STATE_1_VariableLogger;        /* '<S9>/STATE_1' */
  void* STATE_1_ZeroCrossingLogger;    /* '<S9>/STATE_1' */
  void* STATE_1_SampleTimeIdx;         /* '<S9>/STATE_1' */
  void* OUTPUT_1_1_Simulator;          /* '<S9>/OUTPUT_1_1' */
  void* OUTPUT_1_1_SimulationData;     /* '<S9>/OUTPUT_1_1' */
  void* OUTPUT_1_1_DiagnosticManager;  /* '<S9>/OUTPUT_1_1' */
  void* OUTPUT_1_1_VariableLogger;     /* '<S9>/OUTPUT_1_1' */
  void* OUTPUT_1_1_ZeroCrossingLogger; /* '<S9>/OUTPUT_1_1' */
  void* OUTPUT_1_1_SampleTimeIdx;      /* '<S9>/OUTPUT_1_1' */
  void* OUTPUT_1_0_Simulator;          /* '<S9>/OUTPUT_1_0' */
  void* OUTPUT_1_0_SimulationData;     /* '<S9>/OUTPUT_1_0' */
  void* OUTPUT_1_0_DiagnosticManager;  /* '<S9>/OUTPUT_1_0' */
  void* OUTPUT_1_0_VariableLogger;     /* '<S9>/OUTPUT_1_0' */
  void* OUTPUT_1_0_ZeroCrossingLogger; /* '<S9>/OUTPUT_1_0' */
  void* OUTPUT_1_0_SampleTimeIdx;      /* '<S9>/OUTPUT_1_0' */
  int_T STATE_1_Modes;                 /* '<S9>/STATE_1' */
  int32_T STATE_1_MASS_MATRIX_PR;      /* '<S9>/STATE_1' */
  boolean_T STATE_1_CallSimulatorOutput;/* '<S9>/STATE_1' */
  boolean_T OUTPUT_1_1_CallSimulatorOutput;/* '<S9>/OUTPUT_1_1' */
  boolean_T OUTPUT_1_0_CallSimulatorOutput;/* '<S9>/OUTPUT_1_0' */
} DW;

/* Continuous states (default storage) */
typedef struct {
  real_T RC1_ManualLithium_Cell_1RC_equi[4];/* '<S9>/STATE_1' */
} X;

/* State derivatives (default storage) */
typedef struct {
  real_T RC1_ManualLithium_Cell_1RC_equi[4];/* '<S9>/STATE_1' */
} XDot;

/* State disabled  */
typedef struct {
  boolean_T RC1_ManualLithium_Cell_1RC_equi[4];/* '<S9>/STATE_1' */
} XDis;

/* Mass Matrix (global) */
typedef struct {
  int_T ir[2];
  int_T jc[5];
  real_T pr[2];
} MassMatrix;

#ifndef ODE14X_INTG
#define ODE14X_INTG

/* ODE14X Integration Data */
typedef struct {
  real_T *x0;
  real_T *f0;
  real_T *x1start;
  real_T *f1;
  real_T *Delta;
  real_T *E;
  real_T *fac;
  real_T *DFDX;
  real_T *W;
  int_T *pivots;
  real_T *xtmp;
  real_T *ztmp;
  real_T *M;
  real_T *M1;
  real_T *Edot;
  real_T *xdot;
  real_T *fminusMxdot;
  boolean_T isFirstStep;
} ODE14X_IntgData;

#endif

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T I_in;                         /* '<Root>/Current' */
} ExtU;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T V_out;                        /* '<Root>/V_out' */
} ExtY;

/* Real-time Model Data Structure */
struct tag_RTM {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  int_T massMatrixType;
  int_T massMatrixNzMax;
  int_T *massMatrixIr;
  int_T *massMatrixJc;
  real_T *massMatrixPr;
  real_T odeX0[4];
  real_T odeF0[4];
  real_T odeX1START[4];
  real_T odeF1[4];
  real_T odeDELTA[4];
  real_T odeE[4*4];
  real_T odeFAC[4];
  real_T odeDFDX[4*4];
  real_T odeW[4*4];
  int_T odePIVOTS[4];
  real_T odeXTMP[4];
  real_T odeZTMP[4];
  real_T odeMASSMATRIX_M[2];
  real_T odeMASSMATRIX_M1[2];
  real_T odeEDOT[4*4];
  real_T odeXDOT[4];
  real_T odeFMXDOT[4];
  ODE14X_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint16_T clockTick0;
    time_T stepSize0;
    uint16_T clockTick1;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Continuous states (default storage) */
extern X rtX;

/* global MassMatrix */
extern MassMatrix rtMassMatrix;

/* Block signals and states (default storage) */
extern DW rtDW;

/* External inputs (root inport signals with default storage) */
extern ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY rtY;

/* Model entry point functions */
extern void RC1_Manual_initialize(void);
extern void RC1_Manual_step(void);

/* Real-time Model object */
extern RT_MODEL *const rtM;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<Root>/From Workspace' : Unused code path elimination
 * Block '<Root>/Scope' : Unused code path elimination
 * Block '<Root>/Scope1' : Unused code path elimination
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'RC1_Manual'
 * '<S1>'   : 'RC1_Manual/Lithium Cell (1RC equivalent circuit model,  single temperature)'
 * '<S2>'   : 'RC1_Manual/PS-Simulink Converter1'
 * '<S3>'   : 'RC1_Manual/PS-Simulink Converter2'
 * '<S4>'   : 'RC1_Manual/Simulink-PS Converter1'
 * '<S5>'   : 'RC1_Manual/Solver Configuration'
 * '<S6>'   : 'RC1_Manual/PS-Simulink Converter1/EVAL_KEY'
 * '<S7>'   : 'RC1_Manual/PS-Simulink Converter2/EVAL_KEY'
 * '<S8>'   : 'RC1_Manual/Simulink-PS Converter1/EVAL_KEY'
 * '<S9>'   : 'RC1_Manual/Solver Configuration/EVAL_KEY'
 */
#endif                                 /* RTW_HEADER_RC1_Manual_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
