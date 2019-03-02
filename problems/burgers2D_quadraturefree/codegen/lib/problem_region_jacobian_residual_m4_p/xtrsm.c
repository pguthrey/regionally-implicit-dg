/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xtrsm.c
 *
 * Code generation for function 'xtrsm'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "problem_region_jacobian_residual_m4_p.h"
#include "xtrsm.h"

/* Function Definitions */
void b_xtrsm(const double A[32400], double B[180])
{
  int k;
  int kAcol;
  int i;
  for (k = 179; k >= 0; k--) {
    kAcol = 180 * k;
    if (B[k] != 0.0) {
      B[k] /= A[k + kAcol];
      for (i = 0; i < k; i++) {
        B[i] -= B[k] * A[i + kAcol];
      }
    }
  }
}

void xtrsm(const double A[32400], double B[180])
{
  int k;
  int kAcol;
  int i;
  for (k = 0; k < 180; k++) {
    kAcol = 180 * k;
    if (B[k] != 0.0) {
      for (i = k + 1; i + 1 < 181; i++) {
        B[i] -= B[k] * A[i + kAcol];
      }
    }
  }
}

/* End of code generation (xtrsm.c) */
