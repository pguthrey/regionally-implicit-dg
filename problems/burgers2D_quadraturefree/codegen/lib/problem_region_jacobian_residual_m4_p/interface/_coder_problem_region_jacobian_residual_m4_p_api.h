/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_problem_region_jacobian_residual_m4_p_api.h
 *
 * Code generation for function '_coder_problem_region_jacobian_residual_m4_p_api'
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
#include "_coder_problem_region_jacobian_residual_m4_p_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void problem_region_jacobian_residual_m4_p(real_T q[180], real_T q_past
  [180], real_T nuv1, real_T nuv2, real_T solution[180], real_T Jacobian[32400],
  real_T residual[180]);
extern void problem_region_jacobian_residual_m4_p_api(const mxArray * const
  prhs[4], int32_T nlhs, const mxArray *plhs[3]);
extern void problem_region_jacobian_residual_m4_p_atexit(void);
extern void problem_region_jacobian_residual_m4_p_initialize(void);
extern void problem_region_jacobian_residual_m4_p_terminate(void);
extern void problem_region_jacobian_residual_m4_p_xil_terminate(void);

#endif

/* End of code generation (_coder_problem_region_jacobian_residual_m4_p_api.h) */
