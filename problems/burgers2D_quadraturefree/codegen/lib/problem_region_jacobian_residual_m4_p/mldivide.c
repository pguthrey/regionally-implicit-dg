/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mldivide.c
 *
 * Code generation for function 'mldivide'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "problem_region_jacobian_residual_m4_p.h"
#include "mldivide.h"
#include "xtrsm.h"
#include "xgetrf.h"

/* Function Definitions */
void mldivide(const double A[32400], double B[180])
{
  static double b_A[32400];
  int ipiv[180];
  int info;
  double temp;
  memcpy(&b_A[0], &A[0], 32400U * sizeof(double));
  xgetrf(b_A, ipiv, &info);
  for (info = 0; info < 179; info++) {
    if (ipiv[info] != info + 1) {
      temp = B[info];
      B[info] = B[ipiv[info] - 1];
      B[ipiv[info] - 1] = temp;
    }
  }

  xtrsm(b_A, B);
  b_xtrsm(b_A, B);
}

/* End of code generation (mldivide.c) */
