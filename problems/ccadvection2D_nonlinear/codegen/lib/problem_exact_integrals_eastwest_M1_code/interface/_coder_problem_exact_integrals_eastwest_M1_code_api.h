/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_problem_exact_integrals_eastwest_M1_code_api.h
 *
 * Code generation for function '_coder_problem_exact_integrals_eastwest_M1_code_api'
 *
 */

#ifndef _CODER_PROBLEM_EXACT_INTEGRALS_EASTWEST_M1_CODE_API_H
#define _CODER_PROBLEM_EXACT_INTEGRALS_EASTWEST_M1_CODE_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_problem_exact_integrals_eastwest_M1_code_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern real_T problem_exact_integrals_eastwest_M1_code(real_T q, real_T q_past,
  real_T nuv1, real_T nuf);
extern void problem_exact_integrals_eastwest_M1_code_api(const mxArray * const
  prhs[4], int32_T nlhs, const mxArray *plhs[1]);
extern void problem_exact_integrals_eastwest_M1_code_atexit(void);
extern void problem_exact_integrals_eastwest_M1_code_initialize(void);
extern void problem_exact_integrals_eastwest_M1_code_terminate(void);
extern void problem_exact_integrals_eastwest_M1_code_xil_terminate(void);

#endif

/* End of code generation (_coder_problem_exact_integrals_eastwest_M1_code_api.h) */
