/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_problem_region_jacobian_residual_m4_p_mex.c
 *
 * Code generation for function '_coder_problem_region_jacobian_residual_m4_p_mex'
 *
 */

/* Include files */
#include "_coder_problem_region_jacobian_residual_m4_p_api.h"
#include "_coder_problem_region_jacobian_residual_m4_p_mex.h"

/* Function Declarations */
static void c_problem_region_jacobian_resid(int32_T nlhs, mxArray *plhs[3],
  int32_T nrhs, const mxArray *prhs[4]);

/* Function Definitions */
static void c_problem_region_jacobian_resid(int32_T nlhs, mxArray *plhs[3],
  int32_T nrhs, const mxArray *prhs[4])
{
  const mxArray *outputs[3];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 4, 4,
                        37, "problem_region_jacobian_residual_m4_p");
  }

  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 37,
                        "problem_region_jacobian_residual_m4_p");
  }

  /* Call the function. */
  problem_region_jacobian_residual_m4_p_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  problem_region_jacobian_residual_m4_p_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(problem_region_jacobian_residual_m4_p_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  problem_region_jacobian_residual_m4_p_initialize();

  /* Dispatch the entry-point. */
  c_problem_region_jacobian_resid(nlhs, plhs, nrhs, prhs);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_problem_region_jacobian_residual_m4_p_mex.c) */
