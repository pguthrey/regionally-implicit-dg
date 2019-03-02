/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_problem_region_Jacobian_residual_M4_P_api.c
 *
 * Code generation for function '_coder_problem_region_Jacobian_residual_M4_P_api'
 *
 */

/* Include files */
#include "tmwtypes.h"
#include "_coder_problem_region_Jacobian_residual_M4_P_api.h"
#include "_coder_problem_region_Jacobian_residual_M4_P_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131467U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "problem_region_Jacobian_residual_M4_P",/* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[30];
static const mxArray *b_emlrt_marshallOut(const real_T u[30]);
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nuv1,
  const char_T *identifier);
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[30];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *q, const
  char_T *identifier))[30];
static const mxArray *emlrt_marshallOut(const real_T u[900]);
static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[30]
{
  real_T (*y)[30];
  y = e_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static const mxArray *b_emlrt_marshallOut(const real_T u[30])
{
  const mxArray *y;
  const mxArray *m1;
  static const int32_T iv2[1] = { 0 };

  static const int32_T iv3[1] = { 30 };

  y = NULL;
  m1 = emlrtCreateNumericArray(1, iv2, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m1, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m1, iv3, 1);
  emlrtAssign(&y, m1);
  return y;
}

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nuv1,
  const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(nuv1), &thisId);
  emlrtDestroyArray(&nuv1);
  return y;
}

static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = f_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[30]
{
  real_T (*ret)[30];
  static const int32_T dims[2] = { 10, 3 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[30])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *q, const
  char_T *identifier))[30]
{
  real_T (*y)[30];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(q), &thisId);
  emlrtDestroyArray(&q);
  return y;
}

static const mxArray *emlrt_marshallOut(const real_T u[900])
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv0[2] = { 0, 0 };

  static const int32_T iv1[2] = { 30, 30 };

  y = NULL;
  m0 = emlrtCreateNumericArray(2, iv0, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m0, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m0, iv1, 2);
  emlrtAssign(&y, m0);
  return y;
}

static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void problem_region_Jacobian_residual_M4_P_api(const mxArray * const prhs[3],
  int32_T nlhs, const mxArray *plhs[2])
{
  real_T (*Jacobian)[900];
  real_T (*residual)[30];
  real_T (*q)[30];
  real_T (*q_past)[30];
  real_T nuv1;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  Jacobian = (real_T (*)[900])mxMalloc(sizeof(real_T [900]));
  residual = (real_T (*)[30])mxMalloc(sizeof(real_T [30]));

  /* Marshall function inputs */
  q = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "q");
  q_past = emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "q_past");
  nuv1 = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "nuv1");

  /* Invoke the target function */
  problem_region_Jacobian_residual_M4_P(*q, *q_past, nuv1, *Jacobian, *residual);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*Jacobian);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(*residual);
  }
}

void problem_region_Jacobian_residual_M4_P_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  problem_region_Jacobian_residual_M4_P_xil_terminate();
}

void problem_region_Jacobian_residual_M4_P_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

void problem_region_Jacobian_residual_M4_P_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (_coder_problem_region_Jacobian_residual_M4_P_api.c) */
