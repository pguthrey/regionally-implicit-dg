/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_problem_region_Jacobian_residual_M4_P_api.h
 *
 * Code generation for function '_coder_problem_region_Jacobian_residual_M4_P_api'
 *
 */

#ifndef _CODER_PROBLEM_REGION_JACOBIAN_RESIDUAL_M4_P_API_H
#define _CODER_PROBLEM_REGION_JACOBIAN_RESIDUAL_M4_P_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_problem_region_Jacobian_residual_M4_P_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void problem_region_Jacobian_residual_M4_P(real_T q[30], real_T q_past[30],
  real_T nuv1, real_T Jacobian[900], real_T residual[30]);
extern void problem_region_Jacobian_residual_M4_P_api(const mxArray * const
  prhs[3], int32_T nlhs, const mxArray *plhs[2]);
extern void problem_region_Jacobian_residual_M4_P_atexit(void);
extern void problem_region_Jacobian_residual_M4_P_initialize(void);
extern void problem_region_Jacobian_residual_M4_P_terminate(void);
extern void problem_region_Jacobian_residual_M4_P_xil_terminate(void);

#endif

/* End of code generation (_coder_problem_region_Jacobian_residual_M4_P_api.h) */
